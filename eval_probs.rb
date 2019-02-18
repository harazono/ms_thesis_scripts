#! /usr/bin/env ruby
require 'csv'
require 'matrix'
include Math


def count_probs(probs, table, tablesize)
  sum = table.flatten.inject(:+)
  index = tablesize - 1
  for i in [*0..index] do
    for j in [*0..index] do
      probs[i][j] = table[i][j].to_f / sum
    end
  end
end

def print_probs(probs, tablesize)
  index = tablesize - 1
  for i in [*0..index] do
    for j in [*0..index] do
      print sprintf("%10f", probs[j][i].to_f)
      print ", " if j != index
    end
    puts
  end
end

def usage
  print "eval_probs.rb <1mer_ocurrence.csv> <2mer_ocurrence.csv>\n"
end
if ARGV.size < 2
  usage
  exit
end

table_oc1_file = ARGV[0]
table_oc2_file = ARGV[1]

table1 = CSV.table(table_oc1_file, headers:false)
table2 = CSV.table(table_oc2_file, headers:false)

probs1 = Array.new(5).map{Array.new(5, 0)}
probs2 = Array.new(25).map{Array.new(25, 0)}

count_probs(probs1, table1, 5)
count_probs(probs2, table2, 25)


probs2byProbs1 = Array.new(25).map{Array.new(25, 0)}

for i in [*0..24] do
  for j in [*0..24] do
    refidx_head = i / 5
    refidx_tail = i % 5
    qryidx_head = j / 5
    qryidx_tail = j % 5
    probs2byProbs1[i][j] = probs1[refidx_head][qryidx_head].to_f * probs1[refidx_tail][qryidx_tail].to_f
  end
end


print_probs(probs1, 5)
puts
print_probs(probs2byProbs1, 25)
puts
print_probs(probs2, 25)
