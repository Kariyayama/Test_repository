=begin
This script makes all combination of domains to each genes and compare
Input file format is domtblout.
marge make_gene_domain_table.rb and make_table_to_hash.rb and make compare function
USAGE: ruby main.rb INPUTFILE1 INPUTFILE2
=end

class DomCom

  def initialize(group1, group2, com1, com2)
    @gmem1 = group1
    @geme2 = group2
    @domcom1 = com1
    @domcom2 = com2
    @G1geneuniq = Hash.new
    @G2geneuniq = Hash.new
    @G1genecons = Hash.new
    @G2genecons = Hash.new
  end

  def compare(filename1, filename2) # compare domcom1 domain conbination and domcom2 couterpart
    
    ofile1 = File.open("#{filename1.split('/')[-1]}.spec.uniq", "w") # domcom1 specific domain combi
    ofilecon = File.open("#{filename1.split('/')[-1]}_#{filename2.split('/')[-1]}.con.uniq", "w") # common domain conbi
    
    i = 0 # loop counter
    
    domcom2flag = Hash.new
    domcom2.each_key do |dom2key|
      domcom2flag.store(dom2key, 0)
    end
        
    domcom1.each_key do |key| 
      # loop count
      i = count(i)
      
      # main part
      if domcom2.fetch(key, nil) != nil then
        ofilecon.print("#{key.join("\t")}\n") 
        domcom2flag[key] = 1
      else
        ofile1.print("#{key.join("\t")}\n")
      end
      
    end
    
    ofilecon.close
    ofile1.close
    
    # output domcom2 specific domain conbination
    ofile2 = File.open("#{filename2.split('/')[-1]}.spec.uniq", "w")
    domcom2uniq.each_key do |key| 
      if domcom2uniq.fetch(key) == 0 then
        ofile2.print("#{key.join("\t")}\n")
      end
    end
    ofile2.close
    print "Done\n"

  end

  private
  def count(i)
    i.count = 1 + i
    if i.count % 100 == 1 then
      print "#{i.count}.."
    end
    return i.count
  end

end
