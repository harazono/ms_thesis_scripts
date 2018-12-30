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

table = CSV.table(ARGV[0], headers:false)
alt = table.map{ |row| row.map{|col| (log10(col) + 2.5) * 100} }

puts <<"EOS"
#ifndef _SCORE_2MER_MOD_3DIG
#define _SCORE_2MER_MOD_3DIG
int score_2mer_mod_3digit[] = {
EOS
for i in 0..alt.size-1
  for j in 0..alt.size-1
    if alt[j][i].finite? then
      printf("%11d", alt[j][i].round - 100)
    else
#      printf("-inf")
      printf("%11d", -1 * (2 ** 30))
    end
    printf(", ") if j != alt.size-1
  end
#  print "\n"
  if i != alt.size-1 then print ",\n" else print "\n" end
end
puts <<"EOS"
};
#endif
EOS
