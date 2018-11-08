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

ref = Bio::FlatFile.open(Bio::FastaFormat, ARGV[1])
samfile = File.open(ARGV[0], "r")
kmer_mtx = Array.new(5 ** kmer_size).map{Array.new(5 ** kmer_size, 0)}


samfile.each do |line|
	next if line[/^@/]
	sam = SAM.new(line)
	next if sam.rname == "*"
	cigar = CIGAR.new(sam.cigar)
	ref_subseq = String.new
	query_subseq = String.new
	ref_subseq = ref.find{|chr| chr.definition == sam.rname}.seq[sam.pos.to_i, sam.seq.to_s.length]
	p ref_subseq[0, 10]
	next
#	p `ps -o rss= -p #{Process.pid}`.to_i
	query_subseq = sam.seq[cigar.cigarOp[0].length, sam.seq.length - 1]
	ref_aligned = String.new
	query_aligned = String.new
	i = 1
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
	

	i = 0
	while ref_aligned[i + kmer_size] != nil do
		ref_idx = bases2coordinate(ref_aligned[i, kmer_size])
		query_idx = bases2coordinate(query_aligned[i, kmer_size])
		kmer_mtx[ref_idx][query_idx] = kmer_mtx[ref_idx][query_idx] + 1
		i = i + 1
	end
end
