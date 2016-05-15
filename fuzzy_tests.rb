# Universidade Federal da Bahia
# Orientadora: Profa. Dra. Tatiane Nogueira
# Author: Nilton Vasques <nilton.vasques{at}openmailbox.org>
# Description: This script search discover.names files in a path and executes the Flexible 
# Organization Method
# 
# SYNOPSIS
# ruby fuzzy-tests.rb --bases-path PATH [OPTIONS]
#
# OPTIONS
#   --bases-path PATH
#       Sets the PATH, where the script will search for discover.names files.
#
#   --verbose 
#       Sets the script to output verbose logs to STDOUT instead PATH/RESULT_FOLDER/output.log
#       Sets the script to output error logs to STDERR instead PATH/RESULT_FOLDER/error.log
#
#   --timestamp TIMESTAMP
#       Disable the auto generated timestamp that are append in RESULTS_PATH. This is useful
#       when you desire override or execute some steps again.
#
#   --no-cluster
#       Disable clustering and classify steps. Useful when you wants generate gplot and 
#       latex reports.
#
#   --clusters cluster1,cluster2,...,clusterN
#       Sets all clusters methods that will be used in tests. By default, fcm, pcm and pfcm 
#       will be used.
#     
#   --pfcm-descriptors extractor1,extracto2,extractor3
#       Set the pfcm descriptors extractor method
#       Extractors can be: [mixed, tipicality, membership, weighted, mixed-degree] 
#
#   --a DOUBLE
#       Set the pfcm a weight 
#
#   --nd INTEGER 
#       Set the number of descriptors. By default 20 is used
#
#   --b DOUBLE
#       Set the pfcm b weight 
#
#   --ask-number-classes
#       Disable the inference algorithm for number of classes, and request it to user.
#
#   --fuzzyfication M
#       Sets the fuzzyfication value. Default: 2.5
#
#   --fuzzyfication-n N
#       Sets the fuzzyfication value. Default: 2.5
#
#   --error E
#       Sets the minimum error. Default: 0.01 
#
#   --weka-path WEKA_PATH
#       Sets the path where is installed weka.jar. Default: /usr/share/java/weka.jar 
#
# EXAMPLES
#
#    ruby fuzzy-tests.rb --bases-path /home/user/fuzzy/bases/ --clusters fcm,pcm,pfcm --verbose --ask-number-classes
#
#    ruby fuzzy-tests.rb --bases-path /home/user/fuzzy/bases/ 
#
#    ruby fuzzy-tests.rb --bases-path /home/user/fuzzy/bases/ --fuzzyfication 2.0 --error 0.02
#
#    ruby fuzzy-tests.rb --clusters pfcm --pfcm-descriptors mixed,tipicality,membership,weighted,mixed-degree --a 1.0 --b 2.0
#
# MINIMUM REQUIREMENTS
#
#  * Linux
#  * Java 
#  * Ruby >= 2.2.3
#  * gnuplot >= 5.0
#  * latex, pdflatex
#  * libsvm
#  * Weka tool
#!/bin/env/ruby ruby

require 'terminal/progressbar'

$start_time = Time.now

$bases_path = "." 
$results_path = "/home/niltonvasques/Documents/ufba/tcc/results/results-c++-#{Time.now.to_s.gsub(/ |:|-/,"")}" 

$bases = { 
  opinosis: "#{$bases_path}/opinosis.in",
  #newsgroup: "#{$bases_path}/20newsgroup.in",
  #hitech: "#{$bases_path}/hitech.in",
  #nsf: "#{$bases_path}/nsf.in",
  #wap: "#{$bases_path}/wap.in",
  #reuters: "#{$bases_path}/reuters.in", 
}

$bases_class = { 
  opinosis:   3,
  newsgroup:  4,
  hitech:     6,
  nsf:        16,
  wap:        4,
  reuters:    43,
}

def arr_to_hash(arr)
  count  = 0
  arr.reduce({}) { |h,i| h[i] = count; count+=1; h}
end

# DEFAULTS assumptions
$weka_path = "/usr/share/java/weka.jar"
$fuzzyfication = 2.5
$fuzzyfication_n = 2.5
$number_descriptors = 20
$pfcm_a = 1.0
$pfcm_b = 4.0
$error = 0.01
$verbose = false
$no_cluster = false
$mixed_descriptors = "" 
$tipicality_descriptors = ""
$membership_descriptors = ""
$weighted_descriptors = ""
$mixed_degree_descriptors = ""
$clusters = [ "pfcm", "fcm", "pcm" ]
$pfcm_descriptors = [ "tipicality" ]
timestamp = Time.now.to_s.gsub(/ |:|-/,"")

# Handling arguments passed to script
unless ARGV.empty?
  arg_hash = arr_to_hash(ARGV)
  $verbose = ARGV.include? "--verbose"
  $no_cluster = ARGV.include? "--no-cluster"
  $clusters = ARGV[arg_hash["--clusters"]+1].split(",") if ARGV.include? "--clusters"
  $pfcm_descriptors = ARGV[arg_hash["--pfcm-descriptors"]+1].split(",") if ARGV.include? "--pfcm-descriptors"
  $fuzzyfication = ARGV[arg_hash["--fuzzyfication"]+1].to_f if ARGV.include? "--fuzzyfication"
  $number_descriptors = ARGV[arg_hash["--nd"]+1].to_i if ARGV.include? "--nd"
  $error = ARGV[arg_hash["--error"]+1].to_f if ARGV.include? "--error"
  $weka_path = ARGV[arg_hash["--weka-path"]+1] if ARGV.include? "--weka-path"
  timestamp = ARGV[arg_hash["--timestamp"]+1] if ARGV.include? "--timestamp"
  $results_path = "/home/niltonvasques/Documents/ufba/tcc/results/results-c++-#{timestamp}" 
  $mixed_descriptors = "-mixed-desc" if  $pfcm_descriptors.include? "mixed"
  $tipicality_descriptors = "-tip-desc" if  $pfcm_descriptors.include? "tipicality"
  $membership_descriptors = "-memb-desc" if  $pfcm_descriptors.include? "membership"
  $weighted_descriptors = "-mixed-weighted-desc" if  $pfcm_descriptors.include? "weighted"
  $mixed_degree_descriptors = "-mixed-degree-desc" if  $pfcm_descriptors.include? "mixed-degree"
  $fuzzyfication_n = ARGV[arg_hash["--fuzzyfication-n"]+1].to_f if ARGV.include? "--fuzzyfication-n"
  $pfcm_a = ARGV[arg_hash["--a"]+1].to_f if ARGV.include? "--a"
  $pfcm_b = ARGV[arg_hash["--b"]+1].to_f if ARGV.include? "--b"

  if ARGV.include? "--bases-path"
    $bases_path = ARGV[arg_hash["--bases-path"]+1] 
    $local_path = $bases_path 
    $results_path = "#{$bases_path}results-pfcm-fcm-auto#{timestamp}" 

    # REGEX Pattern to try discover how many classes are in a discover.data file.
    INF_CLASS = "cut -f1 -d, | sed -s 's/\\(.*\\/\\).*/\\1/g' | sort | uniq | wc -l"
    $bases = {}
    $bases_class = {}
    unless ARGV.include? "--ask-number-classes"
      print "Searching number of class for each base"
    end
    Dir.glob("#{$bases_path}**/discover.names").map {|i| i.gsub("/discover.names", "")}.each do |i|
      base = i.gsub(/ |\//, "_").gsub(".","")
      $bases[base.to_sym] = i
      unless ARGV.include? "--ask-number-classes"
        print "."
        $bases_class[base.to_sym] = [%x[cat #{i}/discover.data | #{INF_CLASS}].to_i, 3].max
      else
        print "Enter number of classes for #{base}: "
        number_of_classes = STDIN.gets.chomp.to_i
        $bases_class[base.to_sym] = number_of_classes 
      end
    end
    puts ""
  end
end
$verbose_cmd = $verbose ? "": " >> #{$results_path}/output.log 2>> #{$results_path}/error.log"
$verbose_err_cmd = $verbose ? "": " 2>> #{$results_path}/error.log"
$max_depth = 5

# The core function of script, reponsible for delegates actions to right applications
def main
  puts "Creating output folder #{$results_path}"
  system "mkdir #{$results_path}"

  bases_inc = 100.0 / $bases.size
  clusters_inc = bases_inc / $clusters.size

  $bases.each do |k,v|
    $clusters.each do |c| 
      result_folder = "#{$results_path}/#{c}-#{k}"
      puts result_folder if $verbose
      puts "mkdir #{result_folder}"
      raise Extractors unless system "mkdir #{result_folder}"

      puts "Computing optimal clusters through #{k} base using #{c} algorithm"
      raise Exception unless exec_java("#{result_folder}/", c, k, v, $fuzzyfication, $error)

      # Replace class attrs that not fit well in weka
      cmd = ("find #{$results_path} -name \\*.arff -exec sed -i 's/@ATTRIBUTE class\tNUMERIC/@ATTRIBUTE classes\tNUMERIC/g' '{}' \\;#{$verbose_cmd}")
      puts cmd
      system cmd

      puts "Running classifiers on clusters #{k} base"
      if c == "hfcm" or c == "hpcm"
        0.upto($max_depth-1) do |level|
          level_folder = "#{result_folder}-level-#{level}"
          puts level_folder if $verbose
          exec_weka(level_folder)
        end
      else
        case c
        when "fcm"
          result_folder = "#{$results_path}/#{c}-#{k}/clusters.soft-fdcl.arff"
          puts result_folder if $verbose
          exec_weka(result_folder)
        when "pcm"
          result_folder = "#{$results_path}/#{c}-#{k}/clusters.soft-fdcl.arff"
          puts result_folder if $verbose
          exec_weka(result_folder)

          result_folder = "#{$results_path}/#{c}-#{k}/clusters.pdcl.arff"
          puts result_folder if $verbose
          exec_weka(result_folder)
        when "pfcm"
          result_folder = "#{$results_path}/#{c}-#{k}/clusters.soft-fdcl.arff"
          puts result_folder if $verbose
          exec_weka(result_folder)

          result_folder = "#{$results_path}/#{c}-#{k}/clusters.mixed-pdcl.arff"
          puts result_folder if $verbose
          exec_weka(result_folder)
        end
      end
    end # $clusters.each do |c| 
  end
  $end_time = Time.now
  $elapsed = $end_time - $start_time
  days = ( $elapsed / (60*60*24) ).to_i
  hours = ( ($elapsed / (60*60)) % 24 ).to_i
  minutes = ( ($elapsed / (60)) % 60 ).to_i
  seconds = ( $elapsed % 60 ).to_i
  puts "Flexible Documents Organization Method finished in #{days} days #{hours} hours #{minutes} minutes #{seconds} seconds"
end

# Execute the Fuzzy Organization Method implemented in java project
# Params:
# +result_folder+:: Output path for results generated by method
# +cluster+:: The cluster method to be executed (fcm, pcm, pfcm, hfcm, hpcm) 
# +base+:: The base name 
# +base_path+:: The base path
# +m_value+:: The fuzzyfication value
# +error_value+:: The error value 
def exec_java(result_folder, cluster, base, base_path, m_value, error_value)
  method = case cluster
           when "fcm" 
             "-f"
           when "pcm" 
             "-p"
           when "pfcm" 
             "-x"
           end
  fuzzy_jar_cmd = ( "./clustering #{method} -o #{result_folder} -m #{m_value} -e #{error_value} -c #{$bases_class[base]} -r 5 -v  -a #{$pfcm_a} -b #{$pfcm_b} -n #{$fuzzyfication_n} #{$verbose_cmd} < #{base_path}" )

  puts fuzzy_jar_cmd if $verbose
  system fuzzy_jar_cmd
end

# Execute Classify clusters data generated by method in WEKA tool
# Params:
# +result_folder+:: Output path for results generated by method
def exec_weka(result_folder)
  weka_cmd = "./run-weka.sh"
  j48_cmd = "java -classpath #{$weka_path} weka.classifiers.trees.J48 -C 0.25 -M 2 -t #{result_folder} > #{result_folder}J48.weka#{$verbose_err_cmd}"
  nb_cmd = "java -classpath #{$weka_path} weka.classifiers.bayes.NaiveBayes -t #{result_folder} > #{result_folder}NB.weka#{$verbose_err_cmd}"
  nbmult_cmd = "java -classpath #{$weka_path} weka.classifiers.bayes.NaiveBayesMultinomial -t #{result_folder} > #{result_folder}NBMult.weka#{$verbose_err_cmd}"
  svm_cmd = "java -classpath #{$weka_path}:/usr/share/java/libsvm.jar weka.classifiers.functions.LibSVM -S 0 -K 1 -D 3 -G 0.0 -R 0.0 -N 0.5 -M 40.0 -C 2.0 -E 0.001 -P 0.1 -seed 1 -t #{result_folder} > #{result_folder}SVM.weka#{$verbose_err_cmd}"
  knn = []
  1.upto(7) do |x|
    knn[x-1] = "java -classpath #{$weka_path} weka.classifiers.lazy.IBk -K #{x} -W 0 -A \"weka.core.neighboursearch.LinearNNSearch -A \\\"weka.core.EuclideanDistance -R first-last\\\"\" -t #{result_folder} > #{result_folder}KNN#{x}.weka#{$verbose_err_cmd}"
  end

  threads = []
  #WEKA ARFF
  threads << Thread.new do
    puts j48_cmd if $verbose
    system j48_cmd
  end

  threads << Thread.new do
    puts nb_cmd if $verbose
    system nb_cmd
  end

  threads << Thread.new do
    puts nbmult_cmd if $verbose
    system nbmult_cmd
  end

  threads << Thread.new do
    puts svm_cmd if $verbose
    system svm_cmd
  end

  knn.each do |cmd|
    threads << Thread.new do
      puts cmd if $verbose
      system cmd
    end
  end
  threads.each{ |thr| thr.join }
end


# Create a latex comparison table using clusters informations 
# Params:
# +file+:: File object to write data 
# +hfcm+:: hfcm wins 
# +hpcm+:: hpcm wins 
def write_table(file, hfcm, hpcm)
  file.write %q[
      \begin{table}[H]
      \caption{Comparisson}
      \centering
      \begin{tabular}{ |c|c| } 
      \hline
  ]
  file.write "#{$clusters[0]} & #{$clusters[1]} \\\\"
  file.write %q[
      \hline
  ]
  file.write "#{hfcm} & #{hpcm} \\\\\n"
  file.write %q[
      \hline
      \end{tabular}
      \end{table}
  ]
end

# Create a latex bases info table
# Params:
# +file+:: File object to write data 
def documents_table(file)
  file.write "\\section{Bases Info}\n\n"
  file.write %q[
      \begin{table}[H]
      \caption{Bases info}
      \centering
      \begin{tabular}{ |c|c|c|c|c| } 
      \hline
  ]
  file.write "Base & Documents & Terms & Class & N-Grams \\\\"
  $bases.each do |b,v|
    b = b.to_s
    file.write %q[
      \hline
    ]
    terms = %x(cat #{v}/discover.names | wc -l)
    docs = %x(cat #{v}/discover.data | wc -l)
    grams = {}
    lines = File.readlines("#{v}/discover.names")
    lines.each do |line|
      unless line.include?("att_class")
        gram = line.split("_").size 
        grams["#{gram}-gram"] = "true"
      end
    end
    puts v if $verbose
    file.write "#{b.gsub(/_/,"-")} & #{docs} & #{terms} & #{$bases_class[b.to_sym]} & #{grams.keys.map{|k| k.to_s+", "}.join} \\\\\n"
    file.write %q[
      \hline
    ]
  end
  file.write %q[
      \hline
      \end{tabular}
      \end{table}
  ]
end

# Finish latex document
# Params:
# +file+:: File object to write data 
def draw_end(file)
  file.write %s(

    \end{multicols}

    \end{document}
  ) 
end

# Extract hierachy useful data from hfcm and hpcm
# Params:
# +file+:: File object to write data 
# +b+:: The base name
def handle_hierarchy(file, b)
  points = {}
  points[$clusters[0]] = 0
  points[$clusters[1]] = 0

  charts = ""
  5.times do |x|
    img = "#{b}-level-#{x}" 
    charts += "\\includegraphics[width=8cm]{#{img}}\n" if File.exist?("#{img}.png")
    dat_file = File.readlines("#{img}.dat")
    inner_points = {}
    inner_points[$clusters[0]] = 0
    inner_points[$clusters[1]] = 0
    clusters_values = dat_file[0].split(" ")
    dat_file.size.times do |j|
      unless dat_file[j].nil?
        values = dat_file[j].split(" ")
        #    puts values
        if values[1].to_f > values[2].to_f
          inner_points[clusters_values[1]] += 1 
        elsif values[1].to_f < values[2].to_f
          inner_points[clusters_values[2]] += 1 
        end
      end
    end
    value1 = clusters_values[1].nil? ? 0 : inner_points[clusters_values[1]]
    value2 = clusters_values[2].nil? ? 0 : inner_points[clusters_values[2]]
    #puts "#{value1} v1 #{value2} v2"
    if value1 > value2 
      points[clusters_values[1]] += 1 
    elsif value1 < value2 
      points[clusters_values[2]] += 1 
    end
    #puts "#{hfcm_level} hfcm_level x #{hpcm_level} hpcm_level"
  end
  write_table(file, points[$clusters[0]], points[$clusters[1]])
  file.write charts
end

# Extract clustering useful data from fcm, pcm and pfcm
# Params:
# +file+:: File object to write data 
# +b+:: The base name
def handle_clustering(file, b)
  charts = ""
  img = "#{b}" 
  puts "#{img} img"
  charts += "\\includegraphics[width=8cm]{#{img}}\n" if File.exist?("#{img}.png")
  #dat_file = File.readlines("#{$results_path}/#{img}.dat")
  file.write charts
end

# Generate a latex report and plot images using gplot
def create_report
  raise Exception unless system "find #{$results_path}/ -name \\*.png | xargs rm -f "

  $bases.keys.map{|i| i.to_s.gsub(/.*_Testes10/,"")}.each do |b|
    puts "GNUPLOT"
    if %x[which gnuplot].empty?
      puts "ERROR: The images cannot be generated because gnuplot was not found!"
      return
    end
    if $clusters[0] == "hfcm" or $clusters[1] == "hfcm"
      5.times do |x|
        raise Exception unless system "ruby plot.rb #{b}-level-#{x} #{$results_path}/#{$verbose_cmd}"
      end
    else
      raise Exception unless system "ruby plot.rb #{b} #{$results_path}/#{$verbose_cmd}"
    end
  end
  puts "find #{$results_path}/ -name \\*.png -size 0 | xargs rm" if $verbose
  raise Exception unless system "find #{$results_path}/ -name \\*.png -size 0 | xargs rm -f"

  File.open("#{$results_path}/results.tex", "a+") do |file|
    puts "Creating report in latex format..."
    raise Exception unless system "cat article_2.tex | sed -s 's/CURRENTDATE/#{Time.now.strftime("%B %Y")}/g' | sed -s 's/ALGORITHMS/#{$clusters.map{|c| c.upcase }.join(", ")}/g' | tee > #{$results_path}/results.tex #{$verbose_err_cmd}"

    #documents_table(file)

    file.write  "\\begin{multicols}{2} \n"
    $bases.each do |b,path|
      b = b.to_s
      file.write "\\section{#{b.capitalize.gsub(/_/,"-")}}\n\n"
      if $clusters[0] == "hfcm" or $clusters[1] == "hfcm"
        handle_hierarchy(file, b)
        file.write "\\newpage\n\n"
      else
        handle_clustering(file, "#{$results_path}/#{b}")
      end
    end
    draw_end(file)
  end
  if %x[which pdflatex].empty?
    puts "ERROR: The report cannot be generated because pdflatex was not found!"
    return
  end
  raise Exception unless system "pdflatex #{$results_path}/results.tex#{$verbose_cmd}"
  puts "Display report..."
  raise Exception unless system "evince results.pdf#{$verbose_cmd}"
end


if File.exists? $weka_path
  main unless $no_cluster
  create_report
else
  puts "The weka tool was not found in #{$weka_path}."
end
