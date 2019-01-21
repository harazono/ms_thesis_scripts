#! /usr/bin/env ruby
require 'csv'
require 'matrix'
include Math

def usage
  print "array_max-min.rb <table.csv>\n"
end
if ARGV.size < 1
  usage
  exit
end

table = CSV.table(ARGV[0], headers:false)
#alt = table.map{|row| row.map{|col| col != 0 ? (log10(col) * 100).round : -1 * 2 ** 10}}
alt = table.flatten
max = alt.max
min = alt.min

print "max = #{max}\nmin = #{min}\n"
