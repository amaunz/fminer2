# re-implementation of fminer in ruby with integrated last-utils
# args: smi, class 
# TODO mode, minfreq

begin
  require File.dirname(__FILE__) +  "/../last.so"
rescue
  puts "last not found!"
  exit false
end

def read_tab_file(f)
  fc = []
  File.open(f).each do |line|
    line.chomp!
    fc << line.split("\t")
  end
  fc
end

class RubyFminer

  def initialize
    @myFminer = Last::Last.new()
  end

  # Fminer/BBRC re-implementation in Ruby
  # @param[String] the SMI-File holding the molecules
  # @param[String] the CLASS-File holding the activities
  # @param[Integer] minimum frequency
  # @param[Boolean] aromatic perception

  def run_fminer(smi_file, class_file, min_freq, arom=true)
    # Adjust settings
    @myFminer.Reset
    @myFminer.SetConsoleOut(false)
    @myFminer.SetMinfreq(min_freq.to_i)
    @myFminer.SetAromatic(arom)

    # Add compounds below. IMPORTANT! DO NOT CHANGE SETTINGS AFTER ADDING COMPOUNDS!
    smi_h = {} 
    smi=read_tab_file(smi_file)
    smi.each do |s|
      @myFminer.AddCompound(s[1].to_s , s[0].to_i)
      smi_h[s[0].to_i] = s[1].to_s
    end

    cla=read_tab_file(class_file)
    cla.each do |c|
      id=c[0].to_i
      activity=c[2].to_f
      @myFminer.AddActivity(activity , id)
    end

    # gather results for every root node in vector instead of immediate output
    result_str = ""
    (0 .. @myFminer.GetNoRootNodes()-1).each do |j|
      result = @myFminer.MineRoot(j)
      result.each do |res|
        result_str+=res+"\n"
      end
    end
    result_str

  end
end


