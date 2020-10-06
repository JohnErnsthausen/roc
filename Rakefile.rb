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

desc "Git HEAD tarball"
task :tar do
  dirname = Dir.getwd
  dirbasename = File.basename(dirname)
  sh "git archive --format=tar.gz -o #{dirbasename}.tar.gz HEAD"
end
