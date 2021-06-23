#!/usr/bin/ruby

require 'rake/clean'

def method_name
  (/`(.*)'/.match(caller.first) ? $1 : 'Unresolved_Method_Name')
end

files = FileList.new('*', '../src/*')
BLD_FILES = files.collect{ |f| matchData = f.match(/^(.*)\.(o|d)$/); matchData[0] if matchData }.compact
OUT_FILES = [] #files.collect{ |f| matchData = f.match(/^e_(.*)\.out$/); matchData[0] if matchData }.compact
TXT_FILES = [] #files.collect{ |f| matchData = f.match(/^e[_|-](.*)\.txt$/); matchData[0] if matchData }.compact
RUN_EXECS = files.collect{ |f| matchData = f.match(/^e_(.*)$/); matchData[0] if matchData && File.extname(f).empty? && f[-1]!= "." }.compact
GENERATED_FILES = BLD_FILES + OUT_FILES + TXT_FILES + RUN_EXECS

CLEAN.include(GENERATED_FILES)
CLOBBER.include("Makefile")

def generated_files
  GENERATED_FILES.flatten.uniq.collect{|f| f if File.file?(f) }.compact
end

desc "Identify generated files"
task :pf do |t|
  generated_files.each{|f| p f }
end

desc "Remove generated files"
task :remove do |t|
  Rake::Cleaner.cleanup_files(generated_files)
end

desc "Compile an example"
task :compile, [:example] do |t,args|
  args.with_defaults(:example => "foo.c")
  puts "Args used: #{args}"
  begin
    makefile = File.read("RBMakefile.mk").gsub(/##MYEXAMPLE##/,args[:example].ext(''))
    fn = "Makefile"
    File.write(fn, makefile)
    sh "make -j 8"
  rescue Exception => e
    puts "ERROR [#{method_name}]: #{e.to_s} (#{e.class})"
  end
end
