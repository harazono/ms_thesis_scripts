#! /usr/bin/env ruby
# vim: set ts=4 sw=4 :
require 'bio'
require 'optparse'
require '/home/harazono/ms_thesis/ms_thesis_scripts/sam_parser'


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
	STDERR.print "kmer_counter_mod.rb -k <kmer_size> <mapped.sam> <reference.fa> \n"
end


if ARGV.size < 4
	usage
	exit
end

def revcomp(str)
	tmp = str.reverse.tr("ACGTacgt","TGCAtgca")
	return tmp
end



params = ARGV.getopts("k:")
kmer_size = params.find{|k, v| k == "k"}[1].to_i
kmer_size = 2 if kmer_size == 0


#STDERR.print "read reference file : #{ARGV[1]}\n"

ref_raw = Bio::FlatFile.open(Bio::FastaFormat, ARGV[1])
samfile = File.open(ARGV[0], "r")
kmer_mtx = Array.new(5 ** (kmer_size)).map{Array.new(5 ** (kmer_size), 0)}
ref_kmer_mtx = Array.new(5 ** (kmer_size), 0)
ref = REF.new(ref_raw)
#STDERR.print "finish reading reference file\n"
#STDERR.print "reference chromosome are\n==========\n"

ref.refseq.each do |k, v|
# STDERR.print "\"#{k}\"\n"
end
#STDERR.print "==========\n"

###
# kmer_counter design
# count up and fill k-mer * k-mer table.
#
#								 <---reference k-mer--->
#				|
#				|
#				|
#	 query k-mer
#				|
#				|
#				|
#
# ex) 2-mer
# kmer_mtx[i][j]
# i : reference k-mer index
# j : query k-mer index
#
#			 0 1 2 ... i ...
#			 A A A A A C C C C C G G G G G T T T T T * * * * *
#			 A C G T * A C G T * A C G T * A C G T * A C G T *
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
#		GA
#		GC
#		GG
#		GT
#		G*
#		TA
#		TC
#		TG
#		TT
#		T*
#
#
#
#
#
###

line_counter = 1# used for only printing progress
samfile.each do |line|
	next if line[/^@/]
	sam = SAM.new(line)
	next if sam.rname == "*"
	next if sam.seq		== "*"
	cigar					= CIGAR.new(sam.cigar)
	sam_flag			= sam.flag.to_i
	
	#use only primary alignment and supplementary alignment and reverse compliment of them.
	#sample has many supplimentary alignment.
	next if sam_flag != 0 && sam_flag != 16 && sam_flag != 2048 && sam_flag != 2064
	
	ref_subseq		= String.new
	ref_entireseq = String.new
	ref_entireseq = ref.refseq[sam.rname]
	unless ref_entireseq
		STDERR.print "sam.rname = \"#{sam.rname}\"\n"
		puts
		exit
	end
# tail_length = 0
# tail_length = cigar.cigarOp[-1].length if cigar.cigarOp[-1].type == "H"
# ref_subseq = ref_entireseq[sam.pos.to_i - 1, sam.seq.to_s.length.to_i + tail_length]
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
#	 STDERR.print "(c, m, i, d, s, h) = (#{c_sum}, #{m_sum}, #{i_sum}, #{d_sum}, #{s_sum}, #{h_sum})\n"

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
# STDERR.print "ref_subseq		 : \"#{ref_subseq}\"\nquery_subseq	 : \"#{query_subseq}\"\n"
	ref_aligned		= String.new
	query_aligned = String.new
	i = 0
	while cigar.cigarOp[i] != nil do
		cutlen = cigar.cigarOp[i].length
		case cigar.cigarOp[i].type
			when "M" then
				ref_aligned		<< ref_subseq.slice!(0..cutlen - 1)
				query_aligned << query_subseq.slice!(0..cutlen - 1)
			when "I" then
				ref_aligned		<< "-" * cutlen
				query_aligned << query_subseq.slice!(0..cutlen - 1)
			when "D" then
				ref_aligned		<< ref_subseq.slice!(0..cutlen - 1)
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
# STDERR.print "ref_aligned		 : \"#{ref_aligned}\"\nquery_aligned	: \"#{query_aligned}\"\n\n"
	i = 0
	if sam_flag & 16 == 16 then
       ref_aligned   = revcomp(ref_aligned)
       query_aligned = revcomp(query_aligned)
	end
	while ref_aligned[i + kmer_size] != nil	 && query_aligned[i + kmer_size] != nil do
		ref_kmer_str	 = ref_aligned[i, kmer_size]
		query_kmer_str = query_aligned[i, kmer_size]

		ref_idx				 = bases2coordinate(ref_kmer_str)
		query_idx			 = bases2coordinate(query_kmer_str)

		kmer_mtx[ref_idx][query_idx] = kmer_mtx[ref_idx][query_idx] + 1
		ref_kmer_mtx[ref_idx] = ref_kmer_mtx[ref_idx] + 1
		if ref_aligned[i] == "-" && query_aligned[i] == "-" then
			STDERR.print "cigar					 : \"#{sam.cigar}\"\nref_aligned		: ...\"#{ref_aligned[i-1..-1]}\"\nquery_aligned	 : ...\"#{query_aligned[i-1..-1]}\"\n"
			exit
		end

		i = i + 1
	end
	STDERR.print "\rproceed #{line_counter} read"
	line_counter = line_counter + 1
end
STDERR.puts

# print at here
#ref_subseq			: "CGACTATTCC"
#query_subseq		: "CGTCTATTCC"
#ref_aligned		: "CGACTATTCC"
#query_aligned	: "CGTCTATTCC"
#
#							 reference kmer
#							 1,				0,			 0,				0,			 0
# qry kmer		 0,				3,			 0,				0,			 0
#							 0,				0,			 1,				0,			 0
#							 1,				0,			 0,				3,			 0
#							 0,				0,			 0,				0,			 0
#
#
for i in 0..(5 ** kmer_size)-1 do
	for j in 0..(5 ** kmer_size)-1 do
		printf("%8f", (kmer_mtx[j][i] / ref_kmer_mtx[j].to_f / (5 ** kmer_size).to_f))
		printf(", ") if j != (5 ** kmer_size)-1
	end
	print "\n"
end
