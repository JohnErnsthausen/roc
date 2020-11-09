#!/usr/bin/ruby

require 'rake/clean'

SWAP_FILES = FileList.new('**/.*.swp').compact

CLEAN.include("build")
CLOBBER.include("")

BUILD_DIR = "build"
directory BUILD_DIR

def run(cmd)
  Dir.chdir(BUILD_DIR) do
    sh cmd
  end
end

desc "cmake build script"
task :cmake => BUILD_DIR do
  run "cmake .."
end

desc "make"
task :make => :cmake do
  run "make -j 8"
end

desc "make install"
task :install => :make do
  run "sudo make install"
end

desc "make test"
task :test => :make do |t|
  run "make test"
end

desc "run test"
task :run => :make do |t|
  run "./run"
end

desc "Identify swap files"
task :swp do
  SWAP_FILES.each{|f| p f }
end

desc "Remove swap files"
task :rm_swp do
  Rake::Cleaner.cleanup_files(SWAP_FILES)
end

desc "Run clang format on all files"
task :clang do
  matchers = %w{cpp hpp cu c h}.collect{|ext| "**/*.#{ext}" }
  FileList.new(matchers).exclude(/^build/).each{|f| sh "clang-format -i #{f}"}
end

desc "Git HEAD tarball"
task :tar do
  dirname = Dir.getwd
  dirbasename = File.basename(dirname)
  sh "git archive --format=tar.gz -o #{dirbasename}.tar.gz HEAD"
end

desc "Reformat Coeffs"
task :coeffs do
  ofn = File.join(File.dirname(__FILE__), "test", "coeff.txt")
  ifn = File.join(File.dirname(__FILE__), "test", "coeffs.txt")
  lines = File.readlines(ifn)

  File.open(ofn, "w") do |output|
    arr = Array.new
    t = 0
    scale = 0
    lines.each do |line|
      if line.include?('t')
        strarr = line.scan(/\d+\.\d+/)
        if strarr.length==2
          t = strarr[0]
          scale = strarr[1]
          puts "#{line.strip},#{t},#{scale}"
          #output.puts "#{line.strip},#{t},#{scale}"
        else
          puts "WARRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRNING"
        end
        next
      end
      if !line.include?('x')
        arr << line.strip
      else
        output.puts "#{arr.join(' ')} #{scale} #{t}" unless arr.empty?
        arr = Array.new
        next
      end
    end
  end
end
