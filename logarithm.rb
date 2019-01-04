#! /usr/bin/env ruby
require 'csv'
require 'matrix'
include Math

def usage
  print "logarithm.rb <table.csv>\n"
end
if ARGV.size < 1
  usage
  exit
end

diff = 353

table = CSV.table(ARGV[0], headers:false)
alt = table.map{ |row| row.map{|col| log10(col) * 100} }
for i in 0..alt.size-1
  for j in 0..alt.size-1
    if alt[j][i].finite? then
      printf("%d ", alt[j][i] + diff)
    else
      printf("%d ", -1 * 2 ** 10 + diff)
    end
    printf(", ") if j != alt.size-1
  end
  if i != alt.size-1 then print "\n" else print "\n" end
end
