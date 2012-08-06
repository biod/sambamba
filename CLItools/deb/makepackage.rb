#!/usr/bin/env ruby

require 'fileutils'
require 'github_api'
require 'highline'

def print_usage
  puts 'usage: ./makepackage.rb <i386|amd64>'
  puts '       ./makepackage.rb upload <filename>'
  puts
  puts '       first version makes a debian package for specified architecture'
  puts '       second version uploads specified debian package to Github downloads'
end

if ARGV.length.zero?
  print_usage
  exit
end

if ARGV[0] == 'upload'
  if ARGV[1].nil?
    print_usage
    exit
  end
  filename = ARGV[1]
  hl = HighLine.new
  login = hl.ask("GitHub login: ")
  pass = hl.ask("GitHub password: ") {|q| q.echo = '*' }
  github = Github.new :login => login, :password => pass
  response = github.repos.downloads.create 'lomereiter', 'sambamba',
                { :name => filename,
                  :size => File.stat(filename).size,
                  :description => '',
                  :content_type => 'application/x-debian-package'
                }
  github.repos.downloads.upload response, filename
  exit
end

$debdir = Dir.pwd

################################################################################

def get_current_tag
  `git describe --tags \`git rev-list --tags --max-count=1\``.strip
end

def get_current_version
  if $current_version.nil?
    # tag is expected to be of form vX.Y.Z
    tag = get_current_tag
    $current_version = tag[1 .. -1]
  end
  $current_version
end

def generate_debian_binary_file
  File.open('DEBIAN/debian-binary', 'w+') do |f|
    f.puts '6.0'
  end
end

def make_sambamba
  Dir.chdir '../..'
  if $arch == 'i386' then
    `make sambamba-gdc-32`
  else
    `make sambamba-gdc`
  end

  Dir.chdir $debdir
  Dir.chdir $rootdir
end

def copy_sambamba_to_usr_bin
  FileUtils.copy '../../build/sambamba', 'usr/bin'
end

def strip_sambamba_executable
  `strip usr/bin/sambamba`
end

def copy_manpages_to_usr_share_man_man1
  Dir.glob('../../man/*.1').each do |fn|
    FileUtils.copy fn, 'usr/share/man/man1'
  end
end

def gzip_manpages
  `gzip --best usr/share/man/man1/*.1`
end

def get_dependencies
  # TODO: get from dpkg-shlibdeps
  "Depends: libc6 (>= 2.3.2), libgcc1 (>= 1:4.1.1)"
end

def compute_installed_size
  sz = `du -bs .`.to_i
  puts "package size: #{sz}"
  (sz / 1024.0).ceil * 1024
end

def generate_control_file
  text = <<CONTROL
Source: sambamba
Package: sambamba
Version: #{get_current_version}
Section: science
Architecture: #{$arch}
#{get_dependencies}
Installed-Size: #{compute_installed_size}
Maintainer: Artem Tarasov <lomereiter@gmail.com>
Description: Sambamba is a tool for working with SAM/BAM file formats which 
    are widely used in bioinformatics.
    It provides viewing, indexing, sorting, and merging utilities 
    all of which use parallel compression and decompression of BAM files, 
    wisely using resources of modern multi-core machines.
CONTROL
   
  File.open('DEBIAN/control', 'w+') do |f|
    f.puts text
  end
end

def build_debian_package
  Dir.chdir '..'
  `dpkg-deb --build #{$rootdir}`
end

################################################################################

$arch = ARGV[0]

$rootdir = "sambamba-#{get_current_version}_#{$arch}"
FileUtils.remove_entry_secure $rootdir if File.exists? $rootdir
Dir.mkdir $rootdir
Dir.chdir $rootdir
Dir.mkdir 'DEBIAN'
FileUtils.mkdir_p 'usr/bin'
FileUtils.mkdir_p 'usr/share/man/man1'
# stay in rootdir

generate_debian_binary_file
make_sambamba
copy_sambamba_to_usr_bin
strip_sambamba_executable
copy_manpages_to_usr_share_man_man1
gzip_manpages
generate_control_file
build_debian_package
