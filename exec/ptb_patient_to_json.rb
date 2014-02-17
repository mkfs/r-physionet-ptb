#!/usr/bin/env ruby
# Utility to generate a JSON version of the files in an MIT PTB patient data
# directory.
# (c) Copyright 2014 mkfs <https://github.com/mkfs>
# Usage:
#   ptb_patient_to_json.rb dir
# Example:
#   wget -r --cut-dirs=2 -np -nH \
#        http://www.physionet.org/physiobank/database/ptbdb/
#   ptb_patient_to_json.rb ptbdb/patient294

require "rubygems"
require 'json/ext'

=begin rdoc
Return a Hash containing the information encoded in the HEA file header line.
=end
def decode_record_header(line)
  fname, num_sig, samp_freq, num_samp, b_time, b_date = line.strip.split(/\s/)
  fname, num_segs = fname.split('/')

  samp_freq, counter_freq =  samp_freq.split('/')
  counter_freq, counter_base = (counter_freq || '').split('(')
  counter_base, jnk = counter_base.split(')') if counter_base

  { :fname => fname, 
    :num_segments => num_segs ? Integer(num_segs) : nil,
    :num_signals => Integer(num_sig),
    :sample_freq => Float(samp_freq || 250),
    :num_samples => num_samp ? Integer(num_samp) : nil,
    :counter_frequency => counter_freq ? Float(counter_freq) : nil,
    :counter_base_value => counter_base ? Float(counter_base) : nil,
    :base_time => b_time,
    :base_date => b_date
  }
end

=begin rdoc
Return an Array of Hash objects containing the information encoded in each
signal line of a header file.
=end
def decode_signals(lines, num_sig)
  index = -1
  lines.map do |line|
    index += 1

    fname, format, adc_gain, adc_res, adc_zero, initial_value, checksum, 
      block_size, descr = line.strip.split /\s/

    # FORMATxSAMPLES_PER_FRAME:SKEW+OFFSET
    format, samples = format.split('x')
    samples, skew = (samples || '').split(':')
    skew, offset = (skew || '').split('+')

    # ADC_GAIN(BASELINE)/UNITS
    adc_gain, baseline = (adc_gain || '').split('(')
    baseline, units = (baseline || '').split(')/')

    { :filename => fname,
      :index => index,
      :format => Integer(format),
      :samples_per_frame => samples ? Integer(samples) : nil,
      :skew => skew ? Integer(skew) : nil,
      :byte_offset => offset ? Integer(offset) : nil,
      :adc_gain => adc_gain ? Float(adc_gain) : nil,
      :baseline => baseline ? Integer(baseline) : nil,
      :units => units,
      :adc_res => adc_res ? Integer(adc_res) : nil,
      :adc_zero => adc_zero ? Integer(adc_zero) : nil,
      :initial_value => initial_value ? Integer(initial_value) : nil,
      :checksum => checksum ? Integer(checksum) : nil,
      :block_size => block_size ? Integer(block_size) : nil,
      :description => descr
    }
  end
end

=begin rdoc
Return a Hash of information strings stored as comments in the header file.
=end
def decode_info_strings(lines)
  section = 'Info'
  lines.inject({}) do |h,line|
    name, val = line[1..-1].strip.split(':')
    if (! val) || (val.empty?)
      section = name.strip
    else
      h[section] ||= {}
      h[section][name.strip] = val.strip
    end
    h
  end
end

=begin rdoc
Return a Hash containing information in the header file (*.hea).
=end
def read_hea(path)
  return {} if ! File.exist?(path)

  arr = File.read(path)
  lines = arr.lines.reject {|line| (line.start_with? '#') or line.strip.empty?}
  comments = arr.lines.select { |line| line.start_with? '#' }

  h = decode_record_header(lines.shift)
  h[:signals] = decode_signals(lines, h[:num_signals])
  h[:info_strings] = decode_info_strings(comments)

  h
end

=begin rdoc
Generate a Hash [String -> Array] of the 16-bit signed data in each file.
=end
def read_ecg_data_files( data_dir, filenames, signals )
  files = filenames.inject({}) { |h,fname|
    h[fname] ||= signals[fname].first[:format]; h
  }

  files.inject({}) { |h, (name, format)|
    path = File.join(data_dir, name)
    # NOTE: formats besides Format16 are not supported
    next if format != 16
    next if ! File.exist?(path)
    h[name] = File.binread(path).unpack('s*')
    h
  }
end

=begin rdoc
Generate a Hash [String -> Array] containing the ECG data for each lead.
=end
def generate_ecg_lead_data(hdr, filenames, signals, data)
  lead_data = {}
  idx = Hash.new(0)
  hdr[:num_samples].times do |i|
    filenames.each do |fname|
      signals[fname].each do |s|
        lead_data[s[:description]] ||= []
        lead_data[s[:description]] << val = data[fname][idx[fname]]
        idx[fname] += 1
      end
    end
  end

  lead_data
end

=begin rdoc
Read ECG data from directory based on info in header.
=end
def read_dir_ecg_data(data_dir, hdr)
  # list of files containing ECG data
  filenames = hdr[:signals].map { |s| s[:filename] }.uniq

  # list of signals (ECG leads) in each file
  signals = filenames.inject({}) { |h, fname|
    h[fname] = hdr[:signals].select { |s| s[:filename] == fname 
                                     }.sort { |a,b| a[:index] <=> b[:index] }
    h
  }

  data = read_ecg_data_files( data_dir, filenames, signals )
  generate_ecg_lead_data(hdr, filenames, signals, data)
end

=begin rdoc
Generate a Hash containing patient info and ECG data.
=end
def ecg_hash(patient, hdr, lead_data)
  group = File.basename(hdr[:signals].first[:filename])
  id = "#{patient}#{group.gsub(/[^[:digit:]]/, '')}".to_i
  age = hdr[:info_strings]['Info']['age']
  age = (age and age != "n/a") ? Integer(age) : -1
  # TODO : timestamp
  
  { :id => id,
    :group => group,
    :patient => patient,
    :age => age,
    :data => lead_data
  }
end

=begin rdoc
Ensure that ID is unique in array.
=end
def dedupe_ids(arr)
  seen = []
  arr.each do |h|
    while seen.include? h[:id]
      h[:id] += 1
    end
    seen << h[:id]
  end
end

# ----------------------------------------------------------------------
if __FILE__ == $0

  if ARGV.count < 1
    $stderr.puts "Usage: #{$0} DIR [...]"
    $stderr.puts "DIR is a PTB patient data directory (e.g. ptbdb/patient293)"
    exit -1
  end

  arr = ARGV.inject([]) do |arr, path|
    patient = File.basename(path).split('patient').last.to_i

    Dir[ File.join(path, '*.hea') ].each do |hea_path|
      hea = read_hea(hea_path)
      next if ! hea

      dir = File.dirname(hea_path)
      ecg = ecg_hash(patient, hea, read_dir_ecg_data(dir, hea))
      arr << ecg if ecg
    end

    arr
  end

  dedupe_ids(arr)

  puts arr.to_json
end
