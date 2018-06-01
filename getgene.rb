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

  BELONGS = 0
  FILENAME = 1
  CONSERVED = 1
  NON_CONSERVED =0

  def initialize(listfile) # make domain array for each gene  
    list_file = File.open(listfile, "r")  # Input list file
    @group = Hash.new                     # Hash: key => greup, value => [[belongs], [filename]]

    # read input file and make hash gene and domain
    list_file.each_line{ |lst|
      if lst.to_s.include?("#") then # exclude '#' row
      else
        list     = lst.split("\s")
        group    = list[0]
        belongs  = list[BELONGS + 1]
        filename = list[FILENAME + 1]
        if group != nil && belongs != nil && filename != nil then 
          @group.store(group, Array.new(2){Array.new})
          @group.fetch(group)[BELONGS].push(belongs)
          @group.fetch(group)[FILENAME].push(filename)
        end
      end
    }
    @group_number = @group.length
    list_file.close
  end

  def readfile(thrshld)
    @gene_hash   = Hash.new  # key:gene value:domain hash
    @domain_hash = Hash.new  # key:domain value:gene hash
    @gene_belongs = Hash.new # Gene affenity

    @group.each_key do |grp|
      @gene_hash.store(grp, Hash.new)
      @domain_hash.store(grp, Hash.new)
      file_names = @group.fetch(grp)[FILENAME]
      file_names.each {|domtblout| # Input domtblout file
        file = File.open(domtblout, "r") 
        puts "start: make hash #{grp}, #{domtblout}"
        gene_nowgene = [nil, Array.new(2){Array.new}]
        file.each_line{|x|
          gene_nowgene = store_domain(x, thrshld, grp, 
                                      gene_nowgene[0], gene_nowgene[1])
        }
        puts "Done: make hash #{grp}, #{domtblout}"
        file.close
      }
    end
    out_read_file
  end

  def make_domain_hash # domain array to domain conbi hash
    @domcom      = Hash.new # key:domcomb value:gene hash
    @gene_domcom = Hash.new # key:gene value:domcomb hash

    # main part
    @gene_hash.each_key do |grp|
      puts "Start: make combi #{grp}"
      @domcom.store(grp, Hash.new)
      @gene_domcom.store(grp, Hash.new)
      @gene_hash[grp].each_key do |gene_key|
        q = @gene_hash[grp].fetch(gene_key).to_a    
        if q.length > 1 then   # exclude one domain gene
          make_combi(gene_key, grp, q)
        end
      end
      puts "Done: make combi #{grp}"
    end

    out_domcom_to_file # output to file
  end

  def compare_domain(group1, group2)
    exit if can_compare?(group1, group2) == 1 # group1 & group2 exist?
    puts "Start: compare domain"

    # main part
    @dom_mat = Hash.new # Hash: key => domain, value => CONSERVED or NON_CONSERVED
    @domain_hash.each_key do |grp1|
      grp2 = flip(grp1, group1, group2)
      @domain_hash[grp1].each_key do |dom|
        if @domain_hash[grp2].fetch(dom, nil) != nil then
          @dom_mat.store(dom, CONSERVED)
        else
          @dom_mat.store(dom, NON_CONSERVED)
        end
      end
    end
    out_domain_compare_file(group1, group2)
    puts "Done: compare domain"
  end

  def compare_combi(group1, group2) # compare domcom1 domain conbination and domcom2 couterpart
    exit if can_compare?(group1, group2) == 1  # group1 & group2 exist?
    @dom_comb_mat = Hash.new # key => 'domain comb', value => CONSERVED or  NON_CONSERVED 

    # main part
    puts "Start: compare domain conbination"
    @domcom.each_key do |grp1|
      grp2 = flip(grp1, group1, group2)
      @domcom[grp1].each_key do |key| 
        if @domcom[grp2].fetch(key, false) then   @dom_comb_mat.store(key, CONSERVED)
        else  @dom_comb_mat.store(key, NON_CONSERVED)
        end
      end
    end
    out_domcom_compare_file(group1, group2)
    puts "Done: compare domain conbination"
  end

  def make_gene_table
    puts "Gene_name\tCons_Dom\tUniq_Dom\tCons_DomComb\tUniq_DomCom\n"
    
  end

  private
  def out_read_file
    @gene_hash.each_key do |group|
      gene_file = File.open("#{group}.dom", "w")
      @gene_hash[group].each_key do |gene|
        gene_file.print("#{gene}\t#{@gene_hash[group].fetch(gene).join("\t")}\n")
      end
      gene_file.close
    end
    puts "Save gene hash"      
  end

  def out_domcom_to_file # output domain conbinations to outfile
    @domcom.each_key do |group|
      comb_file = File.open("#{group}.comb", "w")
      @domcom[group].each_key do |domain_comb|
        comb_file.print("#{domain_comb.join("\t")}\n")
      end
      comb_file.close
    end
    puts "Save domain conbination hash"
  end

  def out_domain_compare_file(group1, group2)     # output domain conbinations to outfile
    out_file1 = File.open("#{group1}_spec.dom", "w") # domcom1 specific domain combi
    out_file2 = File.open("#{group2}_spec.dom", "w") # domcom2 specific domain combi
    out_file_cnsv = File.open("#{group1}_#{group2}_conserved.dom", "w") # common domain conbi    
    @dom_mat.each_key do |key|
      if @dom_mat[key] == CONSERVED then  
        out_file_cnsv.print("#{key}\n") 
      elsif @dom_mat[key] == NON_CONSERVED && @domain_hash[group1].fetch(key, nil) != nil then 
        out_file1.print("#{key}\n")
      elsif @dom_mat[key] == NON_CONSERVED && @domain_hash[group2].fetch(key, nil) != nil then 
        out_file2.print("#{key}\n")
      else  puts "Unexpected Error in out_domcom_compare_file."; exit
      end
    end
    out_file_cnsv.close
    out_file1.close
    out_file2.close
  end

  def out_domcom_compare_file(group1, group2)     # output domain conbinations to outfile
    out_file1 = File.open("#{group1}_spec.domcom", "w") # domcom1 specific domain combi
    out_file2 = File.open("#{group2}_spec.domcom", "w") # domcom2 specific domain combi
    out_file_cnsv = File.open("#{group1}_#{group2}_conserved.domcom", "w") # common domain conbi    
    @dom_comb_mat.each_key do |key|
      if @dom_comb_mat[key] == CONSERVED then  
        out_file_cnsv.print("#{key.join("\t")}\n") 
      elsif @dom_comb_mat[key] == NON_CONSERVED && @domcom[group1].fetch(key, nil) != nil then
        out_file1.print("#{key.join("\t")}\n")
      elsif @dom_comb_mat[key] == NON_CONSERVED && @domcom[group2].fetch(key, nil) != nil then 
        out_file2.print("#{key.join("\t")}\n")
      else  puts "Unexpected Error in out_domcom_compare_file."; exit
      end
    end
    out_file_cnsv.close
    out_file1.close
    out_file2.close
  end

  def count(i) # counter
    i_count = i + 1
    if i_count % 100 == 1 then print "#{i_count}.."
    end
    return i_count
  end

  def store_domain(line, threshold, grp, gene, nowgene)
    if line.to_s.include?("#") then
    else
      row = line.split("\s")
      pfamid   = row[1]
      geneid   = row[3]
      eval     = row[6]
      alistart = row[17]
      if eval.to_f < threshold.to_f then # threshold E-value
        if @domain_hash[grp].fetch(pfamid, nil) == nil then
          @domain_hash[grp].store(pfamid, Array.new).push(geneid)
        else
          @domain_hash[grp].fetch(pfamid).push(geneid)
        end

        if gene != geneid then # other gene
          # store last gene domain data
          if gene != nil then
            doms = Array.new
            nowgene[1].sort.each do |i|
              doms.push(nowgene[0][nowgene[1].index(i)].split('.')[0])
            end
            @gene_hash[grp].store(gene, doms)
            @gene_belongs.store(gene, [grp, @group.fetch(grp)[BELONGS]])
          end
          # new gene domain memory
          nowgene = Array.new(2){Array.new}
          gene = geneid
          nowgene[0] = pfamid.split()  # pfam accession
          nowgene[1] = alistart.split() # query alignment from

        elsif gene == geneid then # same gene
          nowgene[0].push(pfamid)  # pfam accession
          nowgene[1].push(alistart) # query alignment from
        end
        
      end
    end
    return [gene, nowgene]
  end

  def make_combi(key, group, query)
    @gene_domcom[group].store(key, Array.new)  # store Hash, key=>gene, value=>domain_conbi
    for i in 0..(query.length - 2) do
      for j in (i+1)..(query.length - 1) do
        @gene_domcom[group].fetch(key).push(query[i], query[j])
        @domcom[group].store([query[i], query[j]], 0)
      end
    end
  end

  def can_compare?(group1, group2)
    if @group.fetch(group1, nil) == nil && @group.fetch(group2, nil) == nil then
      puts 'Error: No such groups'
      return 1
    else
      return 0
    end
  end
  
  def flip(now_group, group1, group2)
    if    now_group == group1 then  return group2
    elsif now_group == group2 then  return group1
    else  puts "Error in flip function"; exit
    end
  end

end
