#! /usr/bin/env ruby
require 'bio'
require 'optparse'

class SAM
	attr_accessor :qname, :flag, :rname, :pos, :mapq, :cigar, :rnext, :pnext, :tlen, :seq, :qual
	
	def initialize(line)
		parse(line)
	end

	def parse(line)
		@qname, @flag, @rname, @pos, @mapq, @cigar, @rnext, @pnext, @tlen, @seq, @qual = line.strip.split
	end
end


class CIGAR #CIGARをパーズする。
	class Op #CIGARの要素。M/X/I/Dの情報と長さを保持
		attr_accessor :type, :length
	end
	
	attr_accessor :cigarOp
	
	def initialize(cigar)
		parse(cigar)
	end
	
	def parse(cigar) #CIGAR文字列を順番に読み込んで、OPの配列を作る
		@cigarOp = Array.new
		tmpstring = String.new
		tmpOp = Op.new
		i = 0
		while cigar[i] != nil do
			if cigar[i] =~ Regexp.new("[^MIDNSHPX]") then
				tmpstring << (cigar[i])
			else
				tmpOp.length = tmpstring.to_i
				tmpOp.type = cigar[i].to_s
				@cigarOp.push(tmpOp)
				tmpstring = ""
				tmpOp = Op.new
			end
			i = i + 1
		end
	end
end

class REF

	attr_accessor :refseq
	
	def initialize(ref_raw)
		parse(ref_raw)
	end
	
	def parse(ref_raw)
		@refseq = Hash.new
		ref_raw.each do |f|
			@refseq[f.definition] = f.seq
#			print "@refseq[#{f.definition}] = \"#{@refseq[f.definition]}\"\n"
		end
	end
end


def base2num(base)
	case base 
		when "A", "a" then
			return 0
		when "C", "c" then
			return 1
		when "G", "g" then
			return 2
		when "T", "t" then
			return 3
		when "*" then
			return 4
		else
			return -1
	end
end

def bases2coordinate(bases)
	i = 0
	tmpstr = String.new
	while bases[i] != nil do
		tmpstr << base2num(bases[i]).to_s
		i = i + 1
	end
	return tmpstr.to_i(5)
end



if ARGV.size < 2
	exit
end
params = ARGV.getopts("k:")
kmer_size = params.find{|k, v| k == "k"}[1].to_i
kmer_size = 3 if kmer_size == 0


print "read reference file : #{ARGV[1]}\n"

ref_raw = Bio::FlatFile.open(Bio::FastaFormat, ARGV[1])
samfile = File.open(ARGV[0], "r")
kmer_mtx = Array.new(5 ** kmer_size).map{Array.new(5 ** kmer_size, 0)}

ref = REF.new(ref_raw)
print "finish reading reference file\n"
print "reference chromosome are\n==========\n"

ref.refseq.each do |k, v|
	print "\"#{k}\"\n"
end
print "==========\n"


samfile.each do |line|
	next if line[/^@/]
	sam = SAM.new(line)
	next if sam.rname == "*"
	next if sam.seq   == "*"
	cigar = CIGAR.new(sam.cigar)
	ref_subseq = String.new
	ref_entireseq = String.new
	ref_entireseq = ref.refseq[sam.rname]
	unless ref_entireseq
		print "sam.rname = \"#{sam.rname}\"\n"
		puts
		break
	end
	
	ref_subseq = ref_entireseq[sam.pos.to_i - 1, sam.seq.to_s.length.to_i - 1]
	unless ref_subseq
		p ref_entireseq 
		puts
		break
	end

	query_subseq = String.new
	if cigar.cigarOp[0].type == "S" then
		query_subseq = sam.seq[cigar.cigarOp[0].length, sam.seq.length - 1]
	else
		query_subseq = sam.seq[0, sam.seq.length - 1]
	end
#	print "query_subseq : #{query_subseq}\nref_subseq   : #{ref_subseq}\n"
	ref_aligned = String.new
	query_aligned = String.new
	i = 0
	while cigar.cigarOp[i] != nil do
		cutlen = cigar.cigarOp[i].length
		case cigar.cigarOp[i].type
			when "M" then
				ref_aligned << ref_subseq.slice!(0..cutlen - 1)
				query_aligned << query_subseq.slice!(0..cutlen - 1)
			when "I" then 
				ref_aligned << "*" * cutlen
				query_aligned << query_subseq.slice!(0..cutlen - 1)
			when "D" then
				ref_aligned << ref_subseq.slice!(0..cutlen - 1)
				query_aligned << "*" * cutlen
			when "S" then

			else
				break
		end
		i = i + 1
	end
#	print "#{sam.qname}, #{sam.rname}\n#{ref_aligned[0, 100]}\n#{query_aligned[0, 100]}\n\n"	
	i = 0
	while ref_aligned[i + kmer_size] != nil do
		ref_idx = bases2coordinate(ref_aligned[i, kmer_size])
		query_idx = bases2coordinate(query_aligned[i, kmer_size])
		kmer_mtx[ref_idx][query_idx] = kmer_mtx[ref_idx][query_idx] + 1
		i = i + 1
	end
end


kmer_mtx.each do |f|
	f.each do |g|
		printf("%4d, ", g)
	end
	puts
end

