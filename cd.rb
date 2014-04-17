#!/usr/bin/env ruby

require 'bio'
require 'trollop'

=begin

Dependency: you will have to install imagemagick before doing a bundle install
On fedora: sudo yum install ImageMagick ImageMagick-devel
	
This program takes a multiple sequence alignment and will create a graph showing sequence divergence between pairwise comparisons of sequences chosen by the user
	
=end

#Setup options hash
opts = Trollop::options do
  opt :input, "Mandatory alignment file, currently only ClustalW format", type: :string, short: '-i'   
  opt :seq1, "Mandatory first sequence in comparison", type: :string, short: '-x'
  opt :seq2, "Mandatory second sequence in comparison", type: :string, short: '-y'
  opt :step, "Optional step size", type: :string, short: '-s'
  opt :window, "Optional window size", type: :string, short: '-w'
  opt :genome_groups, "Optional clade grouping file in YAML format.", type: :string, short: '-g'
end


#Check existence of required parameters
abort("Must have an alignment file defined '-h' or '--help' for usage") if opts[:input].nil?
abort("The alignment file must actually exist '-h' or '--help' for usage") unless File.exist?(opts[:input])
abort("Must have an initial sequence defined '-h' or '--help' for usage") if opts[:seq1].nil?
abort("Must have a second sequence defined '-h' or '--help' for usage") if opts[:seq2].nil?

#Update parameters if user set
opts[:step].nil?   ? step = 1      : step = opts[:step].to_i
opts[:window].nil? ? window = 1000 : window = opts[:window].to_i


#Read in alignment
aln = Bio::ClustalW::Report.new(File.read(opts[:input])).alignment
puts aln.class
puts 'Step size set to ' + step.to_s
puts 'Window size set to ' + window.to_s


#Calculate divergence over window between two sequences
def calc_pair_diverg(s1, s2, window, st)
	abort("The size of the two selected sequences differ") if s1.size != s2.size
	scores = []
	puts 'window is ' + window.to_s
	puts 'step is ' + st.to_s
	puts 's1 size is ' + s1.size.to_s
	(window..s1.size).step(st) do |n|
		curr_win_score = 0
		win1 = s1[n-window, window]
		win2 = s2[n-window, window]
		#abort(win1) if win1.size > 1000
		#next
		(0..win1.size).each do |j|
			curr_win_score += 1 if win1[j] == win2[j]
		end
		scores.push(curr_win_score.to_f/window)
	end
	scores
end

scores = calc_pair_diverg(aln[opts[:seq1]],aln[opts[:seq2]], window, step )

scores.each do |num|
    puts num
end