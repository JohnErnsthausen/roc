#!/home/vagrant/.rvm/rubies/ruby-2.2.0/bin/ruby

require 'rake/clean'
require 'csv'

def method_name
  (/`(.*)'/.match(caller.first) ? $1 : 'Unresolved_Method_Name')
end

# $ rake clean
CLEAN.include("*.aux","*.bbl","*.blg","*.log", "*.out", "*.snm", "*.toc", "*.nav", "*.vrb")
# $ rake clobber
CLOBBER.include("*.dvi","*.pdf","images/*.pdf")

desc "Build LaTeX"
task :build, [:texfile] do |t,args|
  args.with_defaults(:texfile => "talk")
  puts "Args used: #{args}"
  tex = "#{args[:texfile]}.tex"
  pdf = "#{args[:texfile]}.pdf"
  flg = "-interaction=batchmode"
  begin
    sh "pdflatex #{flg} #{tex} && pdflatex #{flg} #{tex}"
  rescue Exception => e
    puts "ERROR [#{method_name}]: #{e.to_s} (#{e.class})"
  end
end


