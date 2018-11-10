#! /usr/bin/env ruby


def usage
	print "test_dataset_simulator.rb reference_to_be_made.fa reads_to_be_made.fa\n"
end

def simread(ref)
	#一定の確率でINDELが入るだけ
	#エラー特性は無い
	read = String.new
	i = 0
	while ref[i] != nil do
		case rand(100)
			when 1 then
				tmp = ref[i]
				while tmp == ref[i] do
					tmp = ["A", "C", "G", "T"].sample
				end
				read << tmp
			when 2 then
				read << ref[i] << ["A", "C", "G", "T"].sample
			else
				read << ref[i]
			end
		i = i + 1
	end
	return read
end


if ARGV.size < 2 then
	usage
	exit
end

ref = Array.new(100){["A", "C", "G", "T"].sample}.join

File.open(ARGV[0], "w") do |f|
	f.puts ">testref"
	f.puts ref
end
File.open(ARGV[1], "w") do |f|
	f.puts ">testread"
	f.puts simread(ref)
end

