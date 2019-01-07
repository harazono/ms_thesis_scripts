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
		else
            STDERR.print "base2num got #{base}\n"
			exit(1)
	end
end

def bases2coordinate(bases)
	i = 0
	tmpstr = String.new
	while bases[i] != nil do
		tmpstr << base2num(bases[i]).to_s
		i = i + 1
	end
	return tmpstr.to_i(4)
end

def usage
	STDERR.print "kmer_counter_referenceonly.rb -k <kmer_size> <reference.fa> \n"
end


if ARGV.size < 2
	usage
	exit
end


params = ARGV.getopts("k:")
kmer_size = params.find{|k, v| k == "k"}[1].to_i
kmer_size = 2 if kmer_size == 0

fastafile = ARGV[0]
ref_raw = Bio::FlatFile.open(Bio::FastaFormat, fastafile)
ref_kmer_mtx = Array.new(4 ** (kmer_size), 0)
ref = REF.new(ref_raw)
STDERR.print "input fasta file = #{fastafile}\n"
STDERR.print "k = #{kmer_size}\n"
STDERR.print "finish reading reference file\n"
STDERR.print "reference chromosome\n==========\n"

ref.refseq.each do |k, v|
    STDERR.print "\"#{k}\"\n"
end
STDERR.print "==========\n"


ref.refseq.map{|k, v|
    STDERR.print "begin  counting #{k}\n"
    i = 0
    while v.to_s[i + kmer_size] != nil do
        ref_kmer_str = v.to_s[i, kmer_size]
		ref_idx   = bases2coordinate(ref_kmer_str)
		ref_kmer_mtx[ref_idx] = ref_kmer_mtx[ref_idx] + 1
		i = i + 1
	end
    STDERR.print "finish counting #{k}\n"
}
print "#{ref_kmer_mtx.join("\n")}\n"
STDERR.print "Done.\n"

