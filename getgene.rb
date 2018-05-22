=begin
This script makes all combination of domains to each genes and compare
Input file format is domtblout.
marge make_gene_domain_table.rb and make_table_to_hash.rb and make compare function
USAGE: ruby main.rb INPUTFILE

INPUTFILE format
Group1  SP1  sp1.domtblout
Group1  SP2  sp2.domtblout
.
.
.
Group2  SP1'  sp1'.domtblout
.
.
.

=end

class Domain

  def initialize(listfile) # make domain array for each gene  
    lfile = File.open(listfile, "r")  # Input list file
    @group = Hash.new
    # read input file and make hash gene and domain
    lfile.each_line{|lst|
      list = lst.split("s")
      group    = list[0]
      member   = list[1]
      filename = list[2]
      if @group.fetch(group, nil) == nil then
        @group.store(group, Array.new(2){Array.new})
      end
      @group.fetch(group)[0].push(member)
      @group.fetch(group)[1].push(filename)
    end
    @grp_number = @group.length
  end

  def readfile(thrshld)
    @ghash = Array.new(@grp_number){Hash.new} # key:gene value:domain hash
    @dhash = Array.new(@grp_number){Hash.new} # key:domain value:gene hash
    @gmem  = Hash.new # Gene affenity
    grpnum = 0
    @gene = nil     # gene and domain memory
    @group.each_key do |grp|
      grp[1].each do |file|
        file = File.open(grp[1], "r") # Input domtblout file
        file.each_line{|x|
          store_domain(x, thrshld, grp, grpnum)
        }
        file.close
      end
      grpnum += 1
    end
  end

  def make_domain_hash # domain array to domain conbi hash
    @domcom = Array.new(@grp_number){Hash.new} # key:domcomb value:gene hash
    @gndomcom = Array.new(@grp_number){Hash.new} # key:gene value:domcomb hash
    # main part
    for p in 0..ghash.length do
      @ghash[p].each_key do |key|
        q = @ghash[p].fetch(key).to_a    
        if q.length > 1 then     # exclude one domain gene
          @gndomcomb.store(key, Array.new)
          for i in 0..(q.length - 2) do
            for j in (i+1)..(q.length - 1) do
              @gndomcomb.fetch(key).push(q[i], q[j])
              if @ghash[p].fetch([q[i], q[j]], nil) != nil then
                @domcom[p].fetch([q[i],q[j]]).push(key)
              else
                @domcom[p].store([q[i],q[j]], Array.new).push(key)
              end
            end
          end
        end
      end
    end
  end

  def compare_domain
    @cnsv = Array.new()
    @uniq = Array.new(@group_number){Array.new()}
    for i in 0..1 do
      @dhash[i].each_key do |dom|
        if @dhash[0**i].fetch(dom, nil) != nil then
          if @cnsv.include?(dom) then
          else
            @cnsv.push(dom1)
          end
        else
          @uniq[i].push(dom1)
        end
      end
    end
  end

  def compare_combi(filename1, filename2) # compare domcom1 domain conbination and domcom2 couterpart
    ofile1 = File.open("#{filename1.split('/')[-1]}.spec.uniq", "w") # domcom1 specific domain combi
    ofilecon = File.open("#{filename1.split('/')[-1]}_#{filename2.split('/')[-1]}.con.uniq", "w") # common domain conbi
    
    i = 0 # loop counter
    
    domcom2flag = Hash.new
    domcom[1].each_key do |dom2key|
      domcom2flag.store(dom2key, 0)
    end
        
    domcom[0].each_key do |key| 
      # loop count
      i = count(i)
      
      # main part
      if domcom[1].fetch(key, nil) != nil then
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

  def make_gene_table
    puts "Gene_name\tCons_Dom\tUniq_Dom\tCons_DomComb\tUniq_DomCom\n"
    
  end

  def out_domcom_to_file
    # output domain conbinations to outfile
    for gp in 0..@domcom.length do
      cofile = File.open("Group #{gp}.com.uniq", "w")
      @domcom[p].each_key do |keyd|
        cofile.print("#{keyd.join("\t")}\n")
      end
      cofile.close
    end
    puts "Done: make domain conbination hash"
  end

  private
  def count(i) # counter
    i.count = 1 + i
    if i.count % 100 == 1 then
      print "#{i.count}.."
    end
    return i.count
  end

  def store_domain(line, thr, grp, grpnum)
    if line.to_s.include?("#") then
    else
      row = line.split("\s")
      pfamid   = row[1]
      geneid   = row[3]
      eval     = row[6]
      alistart = row[17]
      if eval.to_f < thr.to_f then # threshold E-value
        if @dhash[grpnum].fetch(pfamid, nil) == nil then
          @dhash[grpnum].store(pfamid, Array.new).push(geneid)
        else
          @dhash[grpnum].fetch(pfamid).push(geneid)
        end

        if @gene != geneid then # other gene
          # store last gene domain data
          if @gene != nil then
            doms = Array.new
            nowgene[1].sort.each do |i|
              doms.push(nowgene[0][nowgene[1].index(i)].split('.')[0])
            end
            @ghash[grpnum].store(gene, doms)
            @gmem.store(gene, [grp, @group.fetch(grp)[0]])
          end
          # new gene domain memory
          nowgene = Array.new(2){Array.new}
          @gene = geneid
          nowgene[0] = pfamid.split()  # pfam accession
          nowgene[1] = alistart.split() # query alignment from
        elsif @gene == geneid then # same gene
          nowgene[0].push(pfamid)  # pfam accession
          nowgene[1].push(alistart) # query alignment from
        end
        
      end
    end
    puts "Done: make hash Group #{grpnum}"   
  end

end
