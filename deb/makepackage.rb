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
  filenames = ARGV[1..-1]
  hl = HighLine.new
  login = hl.ask("GitHub login: ")
  pass = hl.ask("GitHub password: ") {|q| q.echo = '*' }
  github = Github.new :login => login, :password => pass
  filenames.each do |filename|
      response = github.repos.downloads.create 'lomereiter', 'sambamba',
                    { :name => filename,
                      :size => File.stat(filename).size,
                      :description => '',
                      :content_type => 'application/x-debian-package'
                    }
      github.repos.downloads.upload response, filename
  end
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

def generate_copyright_file
  gpl2_notice =<<GPL
 This program is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public
 License as published by the Free Software Foundation; either
 version 2 of the License, or (at your option) any later
 version.
 .
 This program is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the GNU General Public License for more
 details.
 .
 You should have received a copy of the GNU General Public
 License along with this package; if not, write to the Free
 Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 Boston, MA  02110-1301 USA
 .
 On Debian systems, the full text of the GNU General Public
 License version 2 can be found in the file
 `/usr/share/common-licenses/GPL-2'.
GPL

  File.open('usr/share/doc/sambamba/copyright', 'w+') do |f|
    f.puts 'Format: http://www.debian.org/doc/packaging-manuals/copyright-format/1.0/'
    f.puts 'Upstream-Name: sambamba'
    f.puts 'Upstream-Contact: Artem Tarasov <lomereiter@gmail.com>'
    f.puts 'Source: https://github.com/lomereiter/sambamba'
    f.puts
    f.puts 'Files: *'
    f.puts 'Copyright: 2012 Artem Tarasov'
    f.puts 'License: GPL-2+'
    f.puts gpl2_notice
  end
end

def copy_changelog_to_usr_share_doc_sambamba
  FileUtils.copy '../DEBIAN_CHANGELOG', 'usr/share/doc/sambamba/changelog'
end

def gzip_changelog
  `gzip --best usr/share/doc/sambamba/changelog`
end

def make_sambamba
  Dir.chdir '../..'
  if $arch == 'i386' then
    `./build_release_gdc.sh 32`
  elsif $arch == 'amd64' then
    `./build_release_gdc.sh 64`
  else
    puts 'architecture must be one of i386, amd64'
    exit 1
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
Priority: optional
Version: #{get_current_version}
Section: science
Architecture: #{$arch}
#{get_dependencies}
Installed-Size: #{compute_installed_size}
Maintainer: Artem Tarasov <lomereiter@gmail.com>
Description: tool for working with data in SAM/BAM file formats
 It provides viewing, indexing, sorting, and merging utilities 
 all of which use parallel compression and decompression of BAM files, 
 wisely using resources of modern multi-core machines.
CONTROL
   
  File.open('DEBIAN/control', 'w+') do |f|
    f.puts text
  end
end

def add_lintian_override_concerning_embedded_zlib
  File.open('usr/share/lintian/overrides/sambamba', 'w+') do |f|
    f.puts '# Phobos links with Zlib statically'
    f.puts 'sambamba binary: embedded-zlib'
  end
end

def build_debian_package
  Dir.chdir '..'
  `fakeroot dpkg-deb --build #{$rootdir}`
end

################################################################################

$arch = ARGV[0]

$rootdir = "sambamba-#{get_current_version}_#{$arch}"
FileUtils.remove_entry_secure $rootdir if File.exists? $rootdir
Dir.mkdir $rootdir
Dir.chdir $rootdir
Dir.mkdir 'DEBIAN'
FileUtils.mkdir_p 'usr/bin'
FileUtils.mkdir_p 'usr/share/doc/sambamba'
FileUtils.mkdir_p 'usr/share/lintian/overrides'
FileUtils.mkdir_p 'usr/share/man/man1'
# stay in rootdir

generate_copyright_file
copy_changelog_to_usr_share_doc_sambamba
gzip_changelog
make_sambamba
copy_sambamba_to_usr_bin
strip_sambamba_executable
copy_manpages_to_usr_share_man_man1
gzip_manpages
generate_control_file
add_lintian_override_concerning_embedded_zlib
build_debian_package
