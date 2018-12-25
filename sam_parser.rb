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
