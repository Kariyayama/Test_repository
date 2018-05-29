=begin
This script makes all combination of domains to each genes and compare
Input file format is domtblout.
marge make_gene_domain_table.rb and make_table_to_hash.rb and make compare function
USAGE: ruby main.rb INPUTFILE

INPUTFILE format
Group1  SP1  sp1.domtblout
Group1  SP2  sp2.domtblout
...
Group2  SP1'  sp1'.domtblout
...

=end

class Domain

  def initialize(listfile) # make domain array for each gene  
    lfile = File.open(listfile, "r")  # Input list file
    @group = Hash.new
    # read input file and make hash gene and domain
    lfile.each_line{ |lst|
      if lst.to_s.include?("#") then # exclude '#' row
      else
        list = lst.split("\s")
        group    = list[0]
        member   = list[1]
        filename = list[2]
        if group != nil && member != nil && filename != nil then 
          @group.store(group, Array.new(2){Array.new})
          @group.fetch(group)[0].push(member)
          @group.fetch(group)[1].push(filename)
        end
      end
    }
    @grp_number = @group.length
    lfile.close
    # p @group
  end

  def readfile(thrshld)
    @ghash = Array.new(@grp_number){Hash.new} # key:gene value:domain hash
    @dhash = Array.new(@grp_number){Hash.new} # key:domain value:gene hash
    @gmem  = Hash.new # Gene affenity
    @gene = nil     # gene and domain memory
    grpnum = 0
    @group.each_key do |grp|
      # p grp 
      files = @group.fetch(grp)[1]
      files.each {|fl|
        file = File.open(fl, "r") # Input domtblout file
        file.each_line{|x|
          store_domain(x, thrshld, grp, grpnum)
        }
        puts "Done: make hash Group #{grpnum + 1}, #{fl}"   
        file.close
      }
      grpnum += 1
    end
  end

  def make_domain_hash # domain array to domain conbi hash
    @domcom = Array.new(@grp_number){Hash.new} # key:domcomb value:gene hash
    @gndomcom = Array.new(@grp_number){Hash.new} # key:gene value:domcomb hash
    # main part
    for gp in 0..(@grp_number - 1) do
      @ghash[gp].each_key do |key|
        q = @ghash[gp].fetch(key).to_a    
        if q.length > 1 then   # exclude one domain gene
          make_combi(key, gp, q)
        end
      end
      puts "Done, make combi Group #{gp + 1}"
    end
    out_domcom_to_file
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

  def compare_combi # compare domcom1 domain conbination and domcom2 couterpart
    ofile1 = File.open("Group1.spec.uniq", "w") # domcom1 specific domain combi
    ofilecon = File.open("Group1_Group2_conserved.uniq", "w") # common domain conbi    
    i = 0 # loop counter
    
    domcom2flag = Hash.new
    @domcom[1].each_key do |dom2key|
      domcom2flag.store(dom2key, 0)
    end
        
    @domcom[0].each_key do |key| 
      # loop count
      i = count(i)
      
      # main part
      if @domcom[1].fetch(key, nil) != nil then
        ofilecon.print("#{key.join("\t")}\n") 
        domcom2flag[key] = 1
      else
        ofile1.print("#{key.join("\t")}\n")
      end

    end
    ofilecon.close
    ofile1.close
    
    # output domcom2 specific domain conbination
    ofile2 = File.open("Group2.spec.uniq", "w")
    domcom2flag.each_key do |key| 
      if domcom2flag.fetch(key) == 0 then
        ofile2.print("#{key.join("\t")}\n")
      end
    end
    ofile2.close
    print "Done\n"
  end

  def make_gene_table
    puts "Gene_name\tCons_Dom\tUniq_Dom\tCons_DomComb\tUniq_DomCom\n"
    
  end

  private
  def out_domcom_to_file
    # output domain conbinations to outfile
    for gp in 0..(@domcom.length - 1) do
      cofile = File.open("Group_#{gp + 1}.com.uniq", "w")
      @domcom[gp].each_key do |keyd|
        cofile.print("#{keyd.join("\t")}\n")
      end
      cofile.close
    end
    puts "Done: save domain conbination hash"
  end

  def count(i) # counter
    i_count = 1 + i
    if i_count % 100 == 1 then
      print "#{i_count}.."
    end
    return i_count
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
            @nowgene[1].sort.each do |i|
              doms.push(@nowgene[0][@nowgene[1].index(i)].split('.')[0])
            end
            # p doms
            @ghash[grpnum].store(@gene, doms)
            @gmem.store(@gene, [grp, @group.fetch(grp)[0]])
          end
          # new gene domain memory
          @nowgene = Array.new(2){Array.new}
          @gene = geneid
          @nowgene[0] = pfamid.split()  # pfam accession
          @nowgene[1] = alistart.split() # query alignment from

        elsif @gene == geneid then # same gene
          @nowgene[0].push(pfamid)  # pfam accession
          @nowgene[1].push(alistart) # query alignment from
        end
        
      end
    end
  end

  def make_combi(key, group, query)
    @gndomcom[group].store(key, Array.new)
    for i in 0..(query.length - 2) do
      for j in (i+1)..(query.length - 1) do
        @gndomcom[group].fetch(key).push(query[i], query[j])
        if @ghash[group].fetch([query[i], query[j]], nil) != nil then
          @domcom[group].fetch([query[i], qeury[j]]).push(key)
        else
          @domcom[group].store([query[i], query[j]], Array.new).push(key)
        end
      end
    end
  end

end
