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
			@refseq[f.first_name] = f.seq
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
		when "-" then
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

def usage
	#STDERR.print "kmer_counter.rb -k <kmer_size> <mapped.sam> <reference.fa> \n"
end


if ARGV.size < 4
	usage
	exit
end


params = ARGV.getopts("k:")
kmer_size = params.find{|k, v| k == "k"}[1].to_i
kmer_size = 2 if kmer_size == 0


#STDERR.print "read reference file : #{ARGV[1]}\n"

ref_raw = Bio::FlatFile.open(Bio::FastaFormat, ARGV[1])
samfile = File.open(ARGV[0], "r")
kmer_mtx = Array.new(5 ** (kmer_size)).map{Array.new(5 ** (kmer_size), 0)}

ref = REF.new(ref_raw)
#STDERR.print "finish reading reference file\n"
#STDERR.print "reference chromosome are\n==========\n"

ref.refseq.each do |k, v|
#	STDERR.print "\"#{k}\"\n"
end
#STDERR.print "==========\n"

###
# kmer_counter design
# count up and fill k-mer * k-mer table.
#
#                <---reference k-mer--->
#       |
#       |
#       |
#  query k-mer
#       |
#       |
#       |
#
# ex) 2-mer
# kmer_mtx[i][j]
# i : reference k-mer index
# j : query k-mer index
#
#      0 1 2 ... i ...
#      A A A A A C C C C C G G G G G T T T T T * * * * *
#      A C G T * A C G T * A C G T * A C G T * A C G T *
# 0 AA
# 1 AC
# 2 AG
# . AT
# . A*
# . CA
# j CC
# . CG
# . CT
# . C*
#   GA
#   GC
#   GG
#   GT
#   G*
#   TA
#   TC
#   TG
#   TT
#   T*
#
###

samfile.each do |line|
	next if line[/^@/]
	sam = SAM.new(line)
	next if sam.rname == "*"
	next if sam.seq   == "*"
	cigar         = CIGAR.new(sam.cigar)
	ref_subseq    = String.new
	ref_entireseq = String.new
	ref_entireseq = ref.refseq[sam.rname]
	unless ref_entireseq
		STDERR.print "sam.rname = \"#{sam.rname}\"\n"
		puts
		exit
	end
#	tail_length = 0
#	tail_length = cigar.cigarOp[-1].length if cigar.cigarOp[-1].type == "H"
#	ref_subseq = ref_entireseq[sam.pos.to_i - 1, sam.seq.to_s.length.to_i + tail_length]
  c_sum = 0
  m_sum = 0
  i_sum = 0
  d_sum = 0
  s_sum = 0
  h_sum = 0
	i = 0
	while cigar.cigarOp[i] != nil do
		cutlen = cigar.cigarOp[i].length
    c_sum = c_sum + cutlen
		case cigar.cigarOp[i].type
			when "M" then
        m_sum = m_sum + cutlen
			when "I" then
        i_sum = i_sum + cutlen
			when "D" then
        d_sum = d_sum + cutlen
			when "S" then
        s_sum = s_sum + cutlen
      when "H" then
        h_sum = h_sum + cutlen
			else
				break
		end
		i = i + 1
	end
#  STDERR.print "(c, m, i, d, s, h) = (#{c_sum}, #{m_sum}, #{i_sum}, #{d_sum}, #{s_sum}, #{h_sum})\n"

  ref_subseq = ref_entireseq[sam.pos.to_i - 1,sam.pos.to_i + c_sum]
	unless ref_subseq
		STDERR.puts ref_entireseq
		STDERR.puts
		break
	end

	query_subseq = String.new
	if cigar.cigarOp[0].type == "S" then
		query_subseq = sam.seq[cigar.cigarOp[0].length, sam.seq.length]
	else
		query_subseq = sam.seq[0, sam.seq.length]
	end
#	STDERR.print "ref_subseq     : \"#{ref_subseq}\"\nquery_subseq   : \"#{query_subseq}\"\n"
	ref_aligned   = String.new
	query_aligned = String.new
	i = 0
	while cigar.cigarOp[i] != nil do
		cutlen = cigar.cigarOp[i].length
		case cigar.cigarOp[i].type
			when "M" then
				ref_aligned   << ref_subseq.slice!(0..cutlen - 1)
				query_aligned << query_subseq.slice!(0..cutlen - 1)
			when "I" then
				ref_aligned   << "-" * cutlen
				query_aligned << query_subseq.slice!(0..cutlen - 1)
			when "D" then
				ref_aligned   << ref_subseq.slice!(0..cutlen - 1)
				query_aligned << "-" * cutlen
			when "S" then
        #STDERR.print "ref_aligned.length = #{ref_aligned.length} query_aligned.length = #{query_aligned.length}\n" if i != 0
			when "H"
				break
			else
				break
		end
		i = i + 1
	end
  if ref_aligned.length != query_aligned.length then
    STDERR.print "ref_aligned.length = #{ref_aligned.length} query_aligned.length = #{query_aligned.length}\n"
    exit
  end
#	STDERR.print "ref_aligned    : \"#{ref_aligned}\"\nquery_aligned  : \"#{query_aligned}\"\n\n"
	i = 0
	while ref_aligned[i + kmer_size] != nil  && query_aligned[i + kmer_size] != nil do
=begin
		ref_previous   = ref_aligned[i, kmer_size]
		query_previous = query_aligned[i, kmer_size]
		previous_str   = ref_previous + query_previous

		ref_nxt   = ref_aligned[i + kmer_size]
		query_nxt = query_aligned[i + kmer_size]
		nxt_str   = ref_nxt + query_nxt
=end
		ref_kmer_str   = ref_aligned[i, kmer_size]
		query_kmer_str = query_aligned[i, kmer_size]

		ref_idx        = bases2coordinate(ref_kmer_str)
		query_idx      = bases2coordinate(query_kmer_str)

		kmer_mtx[ref_idx][query_idx] = kmer_mtx[ref_idx][query_idx] + 1
    if ref_aligned[i] == "-" && query_aligned[i] == "-" then
      STDERR.print "cigar          : \"#{sam.cigar}\"\nref_aligned    : ...\"#{ref_aligned[i-1..-1]}\"\nquery_aligned  : ...\"#{query_aligned[i-1..-1]}\"\n"
      exit
    end
		i = i + 1
	end
  STDERR.print "finish proccesing #{sam.qname}\n"
end


# print at here
#ref_subseq     : "CGACTATTCC"
#query_subseq   : "CGTCTATTCC"
#ref_aligned    : "CGACTATTCC"
#query_aligned  : "CGTCTATTCC"
#
#              reference kmer
#              1,       0,       0,       0,       0
# qry kmer     0,       3,       0,       0,       0
#              0,       0,       1,       0,       0
#              1,       0,       0,       3,       0
#              0,       0,       0,       0,       0
#
#
for i in 0..(5 ** kmer_size)-1 do
	for j in 0..(5 ** kmer_size)-1 do
		printf("%7d", kmer_mtx[j][i])
		printf(", ") if j != (5 ** kmer_size)-1
	end
	print "\n"
end
=begin
kmer_mtx.each do |f|

	f.each do |g|
		printf("%5d, ", g)
	end
	puts
end
=end
