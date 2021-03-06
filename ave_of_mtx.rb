#! /usr/bin/env ruby
require 'csv'
require 'matrix'
include Math

def usage
  print "ave_of_mtx.rb <table.csv>\n"
end
if ARGV.size < 1
  usage
  exit
end

table = CSV.table(ARGV[0], headers:false)
sum = 0
table.map{ |row| row.map{|col| sum = sum + col } }
p sum
p table.size
p sum.to_f / table.size



=begin
alt = table.map{ |row| row.map{|col| log10(col) + 3} }
for i in 0..alt.size-1
  for j in 0..alt.size-1
    if alt[j][i].finite? then
      printf("%2d", alt[j][i].round)
    else
      printf("%2d", -5)
    end
    printf(", ") if j != alt.size-1
  end
  print "\n"
end
=end
