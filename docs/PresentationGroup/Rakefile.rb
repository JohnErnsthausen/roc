#!!/usr/bin/env rvm 2.7.2

require 'rake/clean'

CLEAN.include("*.aux","*.bbl","*.blg","*.log", "*.out", "*.snm", "*.toc", "*.nav", "*.vrb")
CLOBBER.include("*.dvi","*.pdf","images/*.pdf")

desc "Build LaTeX"
task :build, [:texfile] do |t,args|
  args.with_defaults(:texfile => "talk")
  puts "Args used: #{args}"
  tex = "#{args[:texfile]}.tex"
  pdf = "#{args[:texfile]}.pdf"
  flg = "-interaction=batchmode"
  sh "pdflatex #{flg} #{tex} && pdflatex #{flg} #{tex}"
end


