#!/usr/bin/ruby1.8

# re-implementation of fminer in ruby with integrated last-utils
# args: smi, class 
# TODO mode, minfreq

require 'last'
require 'lu.rb' # AM: need last-utils from http://github.com/amaunz/last-utils!

def read_tab_file(f)
  fc = []
  File.open(f).each do |line|
    line.chomp!
    fc << line.split("\t")
  end
  fc
end


status=false
if $*.size!=4
  status=true
end

if status
  puts "Usage: #{$0} {msa|nls|nop} /path/to/smi-file /path/to/class-file minfreq"
  exit
end

MyFminer = Last::Last.new()
# Adjust settings
MyFminer.SetMaxHops(50)
MyFminer.SetConsoleOut(false)
MyFminer.SetMinfreq($*[3].to_i)


# Add compounds below. IMPORTANT! DO NOT CHANGE SETTINGS AFTER ADDING COMPOUNDS!
smi_h = {} 
smi=read_tab_file($*[1])
smi.each do |s|
  MyFminer.AddCompound(s[1].to_s , s[0].to_i)
  smi_h[s[0].to_i] = s[1].to_s
end

cla=read_tab_file($*[2])
all_hash = Hash.new
cla.each do |c|
  id=c[0].to_i
  activity=c[2].to_f
  MyFminer.AddActivity(activity , id)
  all_hash[id]=activity
end

# gather results for every root node in vector instead of immediate output
xml = ""
(0 .. MyFminer.GetNoRootNodes()-1).each do |j|
  result = MyFminer.MineRoot(j)
  result.each do |res|
    xml << res
  end
end

lu = LU.new 
dom=lu.read(xml) 
smarts=lu.smarts_rb(dom,$*[0]) 
instances=lu.match_rb_hash(smi_h,smarts) 

instances.each do |smarts, ids|
  feat_hash = Hash[*(all_hash.select { |k,v| ids.include?(k) }.flatten)]
  MyFminer.GetRegression() ? p_value = MyFminer.KSTest(all_hash.values, feat_hash.values).to_f : p_value = MyFminer.ChisqTest(all_hash.values, feat_hash.values).to_f 
  print "#{smarts}\t#{p_value.abs}\t#{((p_value>0)?'activating':'deactivating')}\t"
  ids.each do |id|
    print "#{id},"
  end
  puts
end
