require "test/unit"
require 'digest/md5'
require 'yaml'

ENV['FMINER_LAZAR'] = '1'
ENV['FMINER_SMARTS'] = '1'
ENV['FMINER_SILENT'] = '1'

$fminer_file=File.expand_path(File.dirname(__FILE__)) + "/fminer.rb"
begin
  require $fminer_file
rescue Exception
  puts File.new(__FILE__).path + ": file '#{$fminer_file}' not found!"
  exit false
end
$myFminer=RubyFminer.new()

# A Test class for Fminer/BBRC
class TestFminer < Test::Unit::TestCase
  def initialize(foo)
    super(foo)
    configure
  end

  # Determines default configuration, such as paths
  def configure
    @smi_file=File.expand_path(File.dirname(__FILE__)) + "/hamster_carcinogenicity.smi"
    @class_file=File.expand_path(File.dirname(__FILE__)) + "/hamster_carcinogenicity.class"
    @md5_yaml_file=File.expand_path(File.dirname(__FILE__)) + "/fminer_md5.yaml"
    @config=nil
    begin
      @config = YAML.load_file(@md5_yaml_file)
    rescue Exception
      puts File.new(__FILE__).path + ": file '#{@md5_yaml_file}' could not be loaded!"
      exit false
    end
  end

  # Tests default Fminer settings 
  def test_ruby_fminer
    output=$myFminer.run_fminer(@smi_file, @class_file, 2)
    actual_md5=Digest::MD5.hexdigest(output)
    expected_md5=@config['default']
    assert_equal(actual_md5, expected_md5)
  end

  # Tests Kekule representation
  def test_ruby_fminer_kekule
    output=$myFminer.run_fminer(@smi_file, @class_file, 2, false)
    actual_md5=Digest::MD5.hexdigest(output)
    expected_md5=@config['kekule']
    assert_equal(actual_md5, expected_md5)
  end

end
