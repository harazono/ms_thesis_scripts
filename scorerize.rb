#! /usr/bin/env ruby
require 'csv'
require 'matrix'
include Math

def usage
  print "scorerize.rb <table.csv>\n"
end
if ARGV.size < 1
  usage
  exit
end
table = CSV.table(ARGV[0], headers:false)
k = 1
for i in 1..6 do
  if 5 ** i == table[0].size
    k = i
  end
end
diff = 0
diff = -30 + 174 - 30 if k == 1
diff = 200 + 34  - 10 if k == 2
diff = 250 + 133 - 30 if k == 3

alt = table.map{ |row| row.map{|col| log10(col) * 100} }
head_string = "#ifndef _SCORE_#{k}MER_MOD_3DIG\n#define _SCORE_#{k}MER_MOD_3DIG\nint score_#{k}mer_mod_3digit[] = {\n"
puts head_string
for i in 0..alt.size-1
  for j in 0..alt.size-1
    if alt[j][i].finite? then
      printf("%5d", alt[j][i].round + diff)
    else
      printf("%5d", -1 * (2 ** 10) + diff)
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
