
#go over the TODO marked ones
#fix issue so that don't have error or identical in an rows and such, try to add those into another row
#error vector - adding errors there instead of in normal 
#list where to get the panther trees required
#what happens if tree is not found and we skip results? 

# example run :
# Rscript UpdateGraftPoint.r treeGraftoutput2.txt specieslist.csv|specieslist.tsv [saveFolder|savelist.out|savelist.txt]
# Rscript UpdateGraftPoint.r uniprotassignedab70perc2.txt uniprotassignedab70perc2spec.tsv

# can have either names of species list or name of files + name of species list
# it is assumed there are no header names, and input files are in tsv (or csv with a '\t' seperator) format

# output is a text file at the end. 
# AN number is given for new graft point (AN numbers are listed in the PANTHER trees)


onRootFlag = TRUE # change with a switch to make false so we output results

# TODO - add ability to say which figures to create plots for, either giving a number for total plots to make or list of numbers or saying all
# TODO - add switch to include results for if treegrafter puts on root

    args = commandArgs(trailingOnly=TRUE)

    args = c('resgraftupdated6speciesunigraftsdataframe.csv','resgraftupdated6speciesunigraftsdataframeSpecList.csv','new6setupdatedgrafts3.txt') # for subbing in args directly without using Rscript
    #args = c('caljatest1.csv','caljatest2.csv','testname.out')

    library(stringr)
    makeplots = FALSE
    options(warn=2)
    options(scipen = 25)
    ##################################################
    ############# verify/prepare inputs ##############
    ##################################################

    checkBlocks = '-b' %in% args
    noblocks = '-nb' %in% args
    if(noblocks){
        blockFile = vector()
        args = args[-which(args == '-nb')]
    } else if(checkBlocks){ # take the argument after -b and load the file, then delete -b and following arg
        w = which(args == '-b')
        blockFile = args[w+1]
        if(!file.exists(blockFile)){
            print('Block file not found, quitting run')
            quit()
        }
        blockFile = read.csv(blockFile,header=F)
        blockFile = blockFile[,1]
        args = args[-c(w,w+1)]
    } else { # load default blockFile if not specified
        blockFile = read.csv('Data/BlockFile.csv',header=F)[,1]

    }

    # Determine the save location ----------------
    if(length(args) == 3){
        saveloc = args[3]
        #if(grepl('/$',args[3])) {
        #   saveloc = args[3]
        #} else {
        #    saveloc = str_c(args[3],'/')
        #}
        tryCatch(
            {
                write.csv(data.frame(),str_c(saveloc))
            },
            error=function(e) {
                message('Save location invalid')
                quit()
            },
            warning=function(w) {
                message('Save location warning')
                print(w)
            }
        )     
    } else if (length(args) < 2) {
        print('invalid number of inputs, try : Rscript UpdateGraftPoint.r treeGraftoutput.txt specieslist.csv [saveloc|savelist.out|savelist.txt] ')
        quit()
    } else{
        saveloc = './updatedTreeGrafterPoint.out'
    }
    print(paste('On completition file will be saved to : ',saveloc))
    #end save location code ------------------------------------

    #make sure the inputs are correct ------------------------
    ### checking species list can be loaded
    if(!file.exists(args[2])){
        print('Species input file does not exist. Quitting run. ')
        quit()
    }
    if(grepl('.csv$',args[2])) {
        inputDF = read.csv(args[2],header=F)
    } else if(grepl('.tsv$',args[2])) {
        inputDF = read.csv(args[2],sep='\t',header=F)
    } else {
            print('Species input file is not a csv or tsv file. Quitting run.')
            quit()
        }

    ### check if species list has names or just species
    namedCols = FALSE
    if(dim(inputDF)[2] == 2)  {
        inputSpecies = inputDF[,2]
        namedCols = TRUE
        } else if(dim(inputDF)[2] == 1) {
        inputSpecies =  inputDF[,1] 
        } else {
        print('Species Input file has incorrect number of columns. Quitting run. ')
        quit()
    }

    ### checking treegrafter inputs
    if(file.exists(args[1])) {
        treeGrafterRes = readLines(args[1])   
        treeGrafterRes = treeGrafterRes[which(treeGrafterRes!="")] 
    } else{
        print('output results file not found. Quitting run.')
        quit()
    }

    ### making sure the tree grafter file has a longer length
    if(length(treeGrafterRes) > length(inputSpecies)){
        print('Length of species input does not match length of treeGrafter input. Quitting run. ')
        quit()
    }

    #making sure species file is empty
    if(length(inputSpecies) == 0){
        print('No input species in file? Quitting run.')
        quit()
    }

    # getting list of AN numbers and pthr trees from the treeGrafter results
    tgr = str_split(treeGrafterRes,' +|\t|,',simplify=T) # dataframe : id, pthr tree, an, pthrtree:SF# , go terms
    tgr = gsub('"','',tgr) # removing quotation marks
    #print(dim(tgr))
    #print(head(tgr))
    pthr = tgr[,2]
    an = tgr[,3]
    sampleName = tgr[,1]
    unqpthr = unique(pthr) #don't want to reload pthr trees multiple times, so check for each unique panther tree
    idNumber = 1:length(an) #used to keep records in order at end

    if(sum(grepl('AN\\d+',an) == T) != length(an)){
        print('AN list in output file contains non AN number(s)')
        print(an[!grep('AN\\d+',an)])
        quit()
    }
    if(sum(grepl('PTHR\\d\\d\\d\\d\\d',pthr) == T) != length(pthr)){
        print('PTHR list in output file contains non PTHR number(s)')
        quit()
    }

    #TODO fix this, it doesn't work correctly
    if(namedCols){
    mnc = match(inputDF[,1],sampleName)
    if(sum(is.na(mnc)) > 0){
        print('Some sample names in species list do not match sample names in tree grafter output list.')
        quit()
    } else{
        inputSpecies = inputSpecies[mnc]
    }
    }

    ##################################################
    ########### Load in Needed Datasets ##############
    ##################################################

    orthosmrca = readRDS('Data/precomputedOrthologsListFixed.rds') #technically this lists which nodes share a mrca that's a duplication node
    # however there's also ht and unk nodes, so this is actually the paralog list, but used to generate the orthologs with htgroup sets

    #child parent has the ncbi Taxon id numbers for the children and the closest ancestors in the panther taxon tree
    # note if they are in the panther tree the parent is the same as the child
    childparent1 = read.csv('Data/NCBIChildParentOnly1.csv')
    childparent2 = read.csv('Data/NCBIChildParentOnly2.csv')
    childparent = rbind(childparent1,childparent2)
 

    #short names is the 5 letter species designations
    shortnames  = read.csv('Data/ncbitaxonshortnames.csv')[,-1] #remove column of column numbers

    #speciesNamesToID gives the taxon id numbers for the species names in the panther taxon tree
    speciesNamesToID = read.csv('Data/pthrNamesToTaxonIDs.csv')[,-1] 
    #speciesNamesToID = read.csv('Data/pthrNamesToTaxonIDsUpdated.csv')[,-1] 
 
    
    # horizontal transfer groups to differeniate between orthologs and xenologs
    htgroups = readRDS('data/horzTranSubTreeGroups&Unk2.rds')
    htgroups2 = readRDS('Data/horzTranSubTreeGroups&UnkBranchInfo2.rds')






    ####################################################
    #check that blockFile names are in speciesNamesToID#
    ####################################################

        options(scipen = 20) # prevents r from using scientific notation numbering
    if(length(blockFile) > 0){
        if(sum(blockFile %in% speciesNamesToID$id | blockFile %in% speciesNamesToID$name) < length(blockFile)){
            print('Some species names in block file not found in speciesNamesToID')
            print(blockFile[!(blockFile %in% speciesNamesToID$name)])
            quit()
        }
        blockFileBKUP = blockFile
        # Update the blockFile (so they are blocked from adding to ones aboved blocked species, not below or on)
        blockFile = vector()
        for(i in 1:length(blockFileBKUP)){
            blockparent = speciesNamesToID[match(blockFileBKUP[i],speciesNamesToID[,1]),3]
            while(!is.na(blockparent)){
                blockparent = speciesNamesToID[match(blockparent,speciesNamesToID[,3]),4]
                blockFile = c(blockFile,blockparent)
            }
        }
        blockFile = unique(na.omit(blockFile))
    }

    ##################################################
    #Convert Input Species to Closest parent taxon ID#
    ##################################################
    print('converting to closest parent taxon id')
    error = rep(FALSE,length(inputSpecies))
    note = rep('',length(inputSpecies)) #note is the note to be added to the output file, 'error' if 1 or id not recognized
    unqspec = unique(inputSpecies)
    pnames = list() #pnames are the ncbitaxonids of all parents up to LUCA for each input species

    m = match(inputSpecies,shortnames[,1]) #match the 5 letter species names to the short names
    luca = childparent[match(shortnames[m[!is.na(m)],2],childparent[,1]),2] == 1
    if(sum(is.na(m)) > 0) {
        print(paste('ERROR: Some species names not recognized, Example : ',inputSpecies[is.na(m)][1],'- Will be labeled as error in results with no graph output'))
    } 
    if(sum(luca) > 0) print(paste('ERROR: LUCA was the closet ancestor found for input species (Example): ',inputSpecies[luca][1], ' - in taxon tree, these species may be viruses, results marked as error'))
    note[luca] = 'LUCA Closest Ancestor Found'
    error[luca] = TRUE
    note[is.na(m)] = 'Species Name Not Found' # ordering is important here as pnames 1 and '' are changed to both be ''. 
    error[is.na(m)] = TRUE

    m = match(unqspec,shortnames[,1]) #match the 5 letter species names to the short names
    luca = childparent[match(shortnames[m[!is.na(m)],2],childparent[,1]),2] == 1
    for(i in 1:length(unqspec)){
        pnames[[unqspec[i]]] = c(childparent[match(shortnames[m[!is.na(m)][i],2],childparent[,1]),2])
        while(pnames[[unqspec[i]]][length(pnames[[unqspec[i]]])] != 1){
            pnames[[unqspec[i]]] = c(pnames[[unqspec[i]]],speciesNamesToID[match(pnames[[unqspec[i]]][length(pnames[[unqspec[i]]])],speciesNamesToID[,3]),4])
        }
        if(length(pnames[[unqspec[i]]]) == 1) { 
            # set to NA if length 1 
            pnames[[unqspec[i]]] = NA
        }
    }
  

    intAncestors = readRDS('data/internalAncestorsToParent.rds')
    cp = read.csv('data/childparent2025.csv')

    fullpnames = list()
    for(i in 1:length(unqspec)){
        fullpnames[[unqspec[i]]] = shortnames[m[!is.na(m)][i],2]
        while(fullpnames[[unqspec[i]]][length(fullpnames[[unqspec[i]]])] != 1){
            fullpnames[[unqspec[i]]] = c(fullpnames[[unqspec[i]]],cp[match(fullpnames[[unqspec[i]]][length(fullpnames[[unqspec[i]]])],cp[,2]),3])
        }
    }





    ##################################################
    ########### Update The Graft Point ###############
    ##################################################

    newGraft = rep(-1,length(pthr)) #list of the new nodes taken for the graft point
    pathways = vector('list',length(pthr)) #list the pathway between the original graft point and the new graft point for each input
    orthos = vector('list',length(pthr)) #the list of uniprot ids for orthologs for each input
    xenos = vector('list',length(pthr)) # list of paralogs
    paras = vector('list',length(pthr)) # list of xenologs
    unkAmbig = vector('list',length(pthr)) # list of homologs that are ambiguous due to unk 
    unktf = rep(FALSE,length(pthr))
    htAmbig = vector('list',length(pthr)) # list of homologs that are ambiguous due to being horizontal transfers
    htAmbigtf = rep(FALSE,length(pthr))
    orthosPaths = vector('list',length(pthr)) #the list of edges leading to orthologs for each input # Deprecated?
    newAN = rep('',length(pthr))
    distFloat = rep('',length(pthr))
    distBranches = rep('',length(pthr))
    downDup = rep(FALSE,length(pthr))
    upDup = rep(FALSE,length(pthr))
    downHT = rep(FALSE,length(pthr))
    upHT = rep(FALSE,length(pthr))
    #childrenvec = rep('',length(pthr))
    dupUpList = vector('list',length(pthr))
    dupDownList = vector('list',length(pthr))
    htDownList = vector('list',length(pthr))
    htUpList = vector('list',length(pthr))
    inferredMove = rep(FALSE,length(pthr))
    implicitHT = rep(FALSE,length(pthr))
    tgRoot = rep(FALSE,length(pthr))
    leafGraft = rep(0,length(pthr))
    nonAncRoot = rep(FALSE,length(pthr))
    matchBlockList = rep(FALSE,length(pthr))
    exceedBlockList = rep(FALSE,length(pthr))
    leafOrtho = rep(0,length(pthr))
    switchBranch = rep(FALSE,length(pthr))
    jumpUpDist = rep(0,length(pthr))
    dist2.0branches = rep(0,length(pthr))
    equalcounts = 0
    #print(unqpthr)
    htgOrthosNumVec = vector('list',length(pthr))
    hnumVec = rep(-1,length(pthr))
    source('treeGrafterUpgradeFunctionsb.r') #contains the read_panther function to read the tree with additional information, updated from the aphylo package to correct several bugs


    print('start loop')

    # error tracking counters
    errorcounting = 0
    placementcounting = 0
    rooterror = 0
    lucaerror = 0
    



    # the set of unqpthr is used to loop so that each tree is only loaded once
    for(i in 1:length(unique(unqpthr))){ 
        
        #if(i %% 10 == 0) print(paste(placementcounting,errorcounting,rooterror,lucaerror))
        if(i %% 5 == 0) print(paste(i,' out of ',length(unique(unqpthr))))
        if(i < 5) print(paste('On PTHR tree : ',unqpthr[i], ' - ',i,' out of ',length(unique(unqpthr))))
        graftAN = an[pthr == unqpthr[i]] #this is AN number of treeGrafter graft point, graftID for taxonID assigned later
        sampNames = sampleName[pthr == unqpthr[i]] #the assigned name from output file
        idNum = idNumber[pthr == unqpthr[i]] #used to keep track of the records to output in order input
        sampSpec = inputSpecies[pthr == unqpthr[i]] #the input species name
        if(file.exists(str_c('treeFiles/',unqpthr[i],'.tree'))){
            tree = read_panther(str_c('treeFiles/',unqpthr[i],'.tree')) #function loaded from source file
        } else {
            print(paste('Panther Tree not found in data files : ',unqpthr[i],' - Skipping results'))
            orthos[[idNum]] = -5
            next
        }
        uniprot = str_split(tree$tree$tip.label,'UniProtKB=',simplify=T)[,2]
        
        leaves = str_split(tree$tree$tip.label,':',simplify=T)
        anlist = leaves[,1] #AN numbers of leaves
        leaves = str_split(leaves[,2],'\\|',simplify=T) 
        specleaves = leaves[,1] #gives 5 letter species name
        specInternal = tree$internal_nodes_annotations$ancestor #gives full species name (although labeled ancestor should name for given node)
        specIDleaves = speciesNamesToID[match(specleaves,speciesNamesToID[,2]),3] #convert short name to ID number
        specIDinternal = speciesNamesToID[match(specInternal,speciesNamesToID[,1]),3] #converts long name to id number
        #note some internal are NA if they are duplications and thus don't have a species
        specIDinternal[is.na(specIDinternal)] = '-1'
        specIDall = c(specIDleaves,specIDinternal)

        numLeaves = length(specleaves)
        rootpos = numLeaves + 1
        nodeType = tree$internal_nodes_annotations$type
        if(nodeType[1] == 'U') unkAmbig[idNum] = TRUE
        countLeaves = 1:numLeaves
        countInternal = 1:length(specInternal)

        #creating secondary tree with length 1 edges, used to determine the distance between any two nodes, which will be initial method of determining closest new graft point if there are duplicate leaves
        tree2 = tree
        tree2$tree$edge.length = rep(1,length(tree$tree$edge.length))
        treeD = dist.nodes(tree$tree) # distance matrix between any two nodes
        tree2D = dist.nodes(tree2$tree)
        
        # this is to let us easily get count of duplications or HT between any two nodes
        treeDup = tree 
        treeDup$tree$edge.length = rep(0,length(tree$tree$edge.length))
        treeDup$tree$edge.length[treeDup$tree$edge[,1] %in% (which(tree$internal_nodes_annotations$type != 'S')+numLeaves)] = 1
        tree2Dup = dist.nodes(treeDup$tree)

        tree$tree$tip.label = uniprot #changed so labels on graphs just show uniprot id
 
    
        for(j in 1:length(graftAN)){
             jpnames = pnames[[sampSpec[j]]] # list of parent ID numbers for all parent taxons of current samples species
             jpnames.implicit = c(jpnames)
             graftID = '0'
             scheck = 0 
            internalAncestors = specIDinternal %in% jpnames
            newGraftID = jpnames[1]
            placementcounting = placementcounting + 1 # for keeping track of how many graft ids have been placed
            #get graft taxon ID from AN number (from graft point), and get it's position in full node list to be used in the distance matrix
            graftNode = -1
            if(graftAN[j] %in% rownames(tree$internal_nodes_annotations)) { #internal node case
                graftID = speciesNamesToID[match(specInternal[match(graftAN[j],rownames(tree$internal_nodes_annotations))],speciesNamesToID[,1]),3] #get the id number for the specific graft point
                #print(paste('a',i,j))
                if(is.na(graftID)){ #for dealing with treegrafter Assigning to duplication nodes
                    graftID = '-1'
                }
                graftIDPos = numLeaves + match(graftAN[j],rownames(tree$internal_nodes_annotations))
            } else if(graftAN[j] %in% anlist) { #leaf node case
                graftIDPos = match(graftAN[j],anlist)
                graftID = speciesNamesToID[match(specleaves[graftIDPos],speciesNamesToID[,2]),3]
            } else { #not in leaf or internal node????
                print(paste('Error, no taxon id found for initial graft point for : ',sampNames[j]))
                newGraft[idNum[j]] = ''
                error[idNum[j]] = TRUE
                note[idNum[j]] = 'No Taxon for Initial Graft'
                orthos[[idNum[j]]] = -4
                next
            }
            currPos = graftIDPos # this is so we can use currPos later if we don't end up updating it when checking up tree
            lookup = FALSE # to keep track if we look up tree


            # checking to see if implicit along starting branch

            #Comparing the treeGrafter graft node with the closest graph node to input species
            if(graftID == newGraftID){ # This is for the case where we don't update the position because initial graft is already on nearest ancestor
                newGraft[idNum[j]] =  graftIDPos
                newAN[idNum[j]] = graftAN[j]
                pathways[[idNum[j]]] = -1
                currLoc = graftIDPos
                graftNode = graftIDPos
                distBranches[idNum[j]] = 0
                distFloat[idNum[j]] = 0
                orthos[[idNum[j]]] = -1
                next
            } else if( is.na(newGraftID) || (onRootFlag && graftIDPos == rootpos)) { # This case is if pnames had an error in the species name input, hence we don't know where to move it
                newGraft[idNum[j]] = -1
                pathways[[idNum[j]]] = -1
                orthosPaths[[idNum[j]]] = -1
                graftNode = graftIDPos
                distBranches[idNum[j]] = 0
                distFloat[idNum[j]] = 0
                newAN[idNum[j]] = graftAN[j]
                note[idNum[j]] = str_c(note[idNum[j]],'; error with determining ancestors')
                if(is.na(newGraftID)) orthos[[idNum[j]]] = -3
                if(onRootFlag && graftIDPos == rootpos){
                    note[idNum[j]] = 'TreeGrafter Placed on Root'
                    tgRoot[idNum[j]] = TRUE
                    orthos[[idNum[j]]] = -2
                }
                next
            } else{ 

                if(graftID != '-1'){
                    # will use if we don't move down
                    check1 = fullpnames[[sampSpec[j]]] 
                    mspec = match(graftID,speciesNamesToID[,3])
                    if(is.na(mspec)) print('ERROR : mspec - shouldnt be na')
                    if(mspec < 259) {
                        check2 = intAncestors[[mspec]]
                        scheck = sum(check2 %in% check1)
                    }
                }

                ###################################################################################
                # Check if position below current position is more closely related to graft species    
                ###################################################################################
                
                


                updatedGraftpos = graftIDPos
                # check below the initial graft point, see if there is a more closely related node species to input species
                if(graftIDPos >= rootpos) {

         
                    graftChildren = get.descendants(tree, graftIDPos)[-1] # subtract first one as it will be input node
                    graftChildren = graftChildren[graftChildren > numLeaves]
                    graftChildrenTaxID = speciesNamesToID[match(specInternal[graftChildren-numLeaves],speciesNamesToID[,1]),3]
                    # check and see if in the list of jpnames
                    gc = graftChildrenTaxID %in% jpnames
                    if(sum(gc) > 0){
                        w = which(jpnames %in% graftChildrenTaxID[gc])
                        mw = min(w)
                        m = which(graftChildrenTaxID == jpnames[mw])
                        if(length(m) == 1) {
                            updatedGraftpos = graftChildren[m]
                        } else if(length(m) > 1) {
                            #if(i < 100) print(paste('multiple graft children-1',i,j))
                            gcd = c()
                            for(mi in 1:length(m)) gcd = c(gcd,tree2D[graftIDPos,graftChildren[m[mi]]])
                            m = m[which(gcd == min(gcd))]
                            if(length(m) > 1){
                                note[idNum[j]] = str_c('Multiple possible closest taxons, one of multiple with same shortest distance chosen;',note[idNum[j]])
                                m = m[1]
                            } else note[idNum[j]] = str_c('Multiple possible closest taxons, one with shortest distance chosen;',note[idNum[j]])
                            updatedGraftpos = graftChildren[m]
                        }
                        # need to check and see if we position it beneath a duplication node, getting list of all parents of the updated graft point
                    }   
                }

                    
                
                
                
                # Perform updates if we updated graft position going down
                if(updatedGraftpos != graftIDPos){
                    # finding number of duplication and ht nodes between initial and updated position
                    finalPos = updatedGraftpos
                    pathBetween = c(updatedGraftpos)
                    while(TRUE){
                        pathBetween = c(pathBetween,tree$tree$edge[tree$tree$edge[,2] == pathBetween[length(pathBetween)],1])
                        if(pathBetween[length(pathBetween)] == graftIDPos) break
                    }
                    typeTemp = nodeType[pathBetween-numLeaves]
                    if(sum(typeTemp == 'D' ) > 0){
                        downDup[idNum[j]] = TRUE
                        dupDownList[[idNum[j]]] = pathBetween[typeTemp == 'D']
                    }
                    if(sum(typeTemp == 'T' ) > 0){
                        downHT[idNum[j]] = TRUE
                        htDownList[[idNum[j]]] = pathBetween[typeTemp == 'T']
                    }
                    # Assigning other statisics 
                    newAN[idNum[j]] = tree$tree$node.label[updatedGraftpos-numLeaves]
                    newGraft[idNum[j]] = updatedGraftpos
                    distBranches[idNum[j]] = tree2D[graftIDPos,updatedGraftpos]
                    distFloat[idNum[j]] = treeD[graftIDPos,updatedGraftpos]                 
                    graftNode = updatedGraftpos
                    #get pathway between the two graft points
                    currLoc = graftNode #note graftNode is position in node list of updated graft point
                    pathway = pathBetween # inclusive of initial and final
                    graftCurr = graftIDPos
                } else if (scheck > 0){
                    print(paste('scheck prevented us from moving in this instance',i,j,scheck,graftID,sampSpec[j]))
                    newGraft[idNum[j]] =  graftIDPos
                    newAN[idNum[j]] = graftAN[j]
                    pathways[[idNum[j]]] = -1
                    currLoc = graftIDPos
                    graftNode = graftIDPos
                    distBranches[idNum[j]] = 0
                    distFloat[idNum[j]] = 0
                    orthos[[idNum[j]]] = -1
                    next
                }

                
                else { # looking up if we didn't assign a graft position
                    
                ###################################################################################
                # Looking up the tree    
                ###################################################################################




                print(paste(i,j))
                #pnames.branch = c(branchimplict[idNum[j]],jpnames)
                    lookup = TRUE
                    parentVec = c(graftIDPos)
                    while(TRUE){  #list of all parents of the original graft point
                        parentVec = c(parentVec,tree$tree$edge[tree$tree$edge[,2] == parentVec[length(parentVec)],1])
                        if(parentVec[length(parentVec)] == rootpos ) break
                    }
                    # remove leaves
                    if(parentVec[1] <= numLeaves) {
                        leafGraft[idNum[j]] = parentVec[1]
                    }
                    parentVec = parentVec[parentVec > numLeaves]
                    
                    # subtract numLeaves
                    parentVec = parentVec - numLeaves
                    # convert the vecs to ids
                    parentVecID = speciesNamesToID[match(specInternal[parentVec],speciesNamesToID[,1]),3] # some are NA due to duplications
                    # then after you get the list of ids for graft point, do a which in the list of children
                    parentType = nodeType[parentVec]

                    w = which(parentVecID %in% jpnames) #.branch)
                    wblock = which(specIDinternal[parentVec] %in% blockFile) 
                    mw = 0   
                    if(length(w) > 0) mw = min(w) # this is the position of closest related parent of graft point that is also in the list of parents of the input species
                    mwblock = Inf
                    if(length(wblock)>0) mwblock = min(wblock)
                    #if(mw == mwblock) { # closest ancestor in block list
                    #    equalcounts = equalcounts + 1 # just to see how many end up going to a top level node that is blocked (so could potentially go further if not for block)
                    #    matchBlockList[idNum[j]] = TRUE
                    #}
                    #if(length(w) == 0 | mw > mwblock) { # the case if there are no valid parents or if we are blocked?
                    sumblock = 0
                    sumblock = sum(blockFile %in% specIDinternal) # don't allow jump up to root if blockfile term in tree (may not be at root so may sneak behind otherwise)
                    if( mw >= mwblock || (sumblock > 0 && length(w) == 0 )) { # if the closest ancestor above us is in a blocking group
                        print('blocked --------------------------------------------')
                        #print(paste('Error, no parent of graft point in list of parents for input species : ',sampNames[j]))
                        exceedBlockList[idNum[j]] = TRUE
                        newGraft[idNum[j]] = -1
                        newAN[idNum[j]] = graftAN[j]
                        error[idNum[j]] = TRUE
                        pathways[[idNum[j]]] = -1
                        orthos[[idNum[j]]] = -8
                        orthosPaths[[idNum[j]]] = -1
                        graftNode = graftIDPos
                        distBranches[idNum[j]] = 0
                        distFloat[idNum[j]] = 0   
                        note[idNum[j]] = str_c(note[idNum[j]],'Block File Prevents Placement, Original Position Retained')   
                        next
                    }

                    
                    if(length(w) == 0 ) {
                        if(sumblock == 0){  
                            currPos = rootpos
                            nonAncRoot[idNum[j]] = TRUE
                            mw = length(parentVec)
                            jumpUpDist[idNum[j]] =  mw
                        }
                    } else {
                        currPos = parentVec[mw] + numLeaves
                        jumpUpDist[idNum[j]] =  mw

                    }
                    if(sumblock == 0){
                        if(sum(parentType[1:mw] == 'D') > 0) {
                            upDup[idNum[j]] = TRUE
                            dupUpList[[idNum[j]]] = parentVec[1:mw][parentType[1:mw] == 'D']
                        }
                        if(sum(parentType[1:mw] == 'T') > 0) {
                            upHT[idNum[j]] = TRUE
                            htUpList[[idNum[j]]] = parentVec[1:mw][parentType[1:mw] == 'T']  
                        }
                    }
                    finalPos = currPos # Final position is updated to final value, currPos stops after moving up being updated

                    ###################################################################################
                    # Looking down the tree once more from updated position
                    ###################################################################################

                    if(currPos >= rootpos) { # this should always be true
                        graftChildren = get.descendants(tree, currPos)[-1] # subtract first one as it will be input node
                        graftChildren = graftChildren[graftChildren > numLeaves]
                        graftChildrenTaxID = speciesNamesToID[match(specInternal[graftChildren-numLeaves],speciesNamesToID[,1]),3]
                        # check and see if in the list of jpnames
                        gc = graftChildrenTaxID %in% jpnames
                        if(sum(gc) > 0){
                            w = which(jpnames %in% graftChildrenTaxID[gc])
                            mw = min(w)
                            m = which(graftChildrenTaxID == jpnames[mw])
                            if(length(m) == 1) {
                                finalPos = graftChildren[m]
                            } else if(length(m) > 1) {
                                #if(i < 100) print(paste('multiple graft children-1',i,j))
                                gcd = c()
                                for(mi in 1:length(m)) gcd = c(gcd,tree2D[graftIDPos,graftChildren[m[mi]]])
                                m = m[which(gcd == min(gcd))]
                                if(length(m) > 1){
                                    note[idNum[j]] = str_c('Multiple possible closest taxons, one of multiple with same shortest distance chosen;',note[idNum[j]])
                                    m = m[1]
                                } else note[idNum[j]] = str_c('Multiple possible closest taxons, one with shortest distance chosen;',note[idNum[j]])
                                finalPos = graftChildren[m]
                            }
                            # need to check and see if we position it beneath a duplication node, getting list of all parents of the updated graft point
                    }
                }
                    
                # Perform updates if we updated graft position going down
                if(finalPos != currPos){
                    # finding number of duplication and ht nodes between initial and updated position
                    pathBetween = c(finalPos)
                    while(TRUE){
                        pathBetween = c(pathBetween,tree$tree$edge[tree$tree$edge[,2] == pathBetween[length(pathBetween)],1])
                        if(pathBetween[length(pathBetween)] == currPos) break
                    }
                    typeTemp = nodeType[pathBetween-numLeaves]
                    if(sum(typeTemp == 'D' ) > 0){
                        downDup[idNum[j]] = TRUE
                        dupDownList[[idNum[j]]] = c(dupDownList[[idNum[j]]],pathBetween[typeTemp == 'D'])
                    }
                    if(sum(typeTemp == 'T' ) > 0){
                        downHT[idNum[j]] = TRUE
                        htDownList[[idNum[j]]] = c( htDownList[[idNum[j]]],pathBetween[typeTemp == 'T'])
                    }
                }

                if(nonAncRoot[idNum[j]] && finalPos == rootpos){
                    finalPos = graftIDPos
                    currPos = graftIDPos
                }
                if(finalPos < rootpos && finalPos != graftIDPos) print('ERROR : final position leaf other than initial graft')
                if(finalPos >= rootpos){
                    # Assigning other statisics 
                    testpos = 1
                    newAN[idNum[j]] = tree$tree$node.label[finalPos-numLeaves]
                    newGraft[idNum[j]] = finalPos
                    distBranches[idNum[j]] = tree2D[graftIDPos,finalPos]
                    distFloat[idNum[j]] = treeD[graftIDPos,finalPos]                 
                    graftNode = finalPos
                    #get pathway between the two graft points
                    currLoc = graftNode #note graftNode is position in node list of updated graft point
                    #pathway = pathBetween # inclusive of initial and final
                    graftCurr = graftIDPos
                    
                    }
                    
                }

                

            }
        
            
           

            ##################################
            #######  Check Positioning #######
            ##################################
            # This is to deal with various cases based on updated graft point


            # DISABLED - this is only the case if treegrafter puts on root and in that case we don't output unless we have flag
            # if we end up on root from moving up move back to original graft point
            # two cases, either root was an ancestor or wasn't
            #caseA = finalPos == rootpos && rootpos != graftIDPos # check if we ended at root
            if(FALSE){ # was if(caseA)
                newAN[idNum[j]] = tree$tree$node.label[graftIDPos-numLeaves]
                newGraft[idNum[j]] = graftIDPos
                distBranches[idNum[j]] = 0
                distFloat[idNum[j]] = 0               
                graftNode = graftIDPos
                orthos[[idNum[j]]] = -9
                #get pathway between the two graft points
                currLoc = graftIDPos #note graftNode is position in node list of updated graft point
                pathway = vector() # inclusive of initial and final
                graftCurr = graftIDPos
                downDup[idNum[j]] = FALSE
                dupDownList[[idNum[j]]] = -1
                downHT[idNum[j]] = FALSE
                htDownList[[idNum[j]]] = -1
                upDup[idNum[j]] = FALSE
                dupUpList[[idNum[j]]] = -1
                upHT[idNum[j]] = FALSE
                htUpList[[idNum[j]]] = -1

           }
            # not sure what this was meant for
            #if(finalPos == rootpos && rootpos != graftIDPos){
            #    note[idNum[j]] = 'TreeGrafter Placed on Root'                
            #}




            # There is an issue here in that we treat edge distance and explicit nodes as mattering more for choosing to go through an ht going down
            # here we are treating implicit as mattering more

            # We also don't want to move through an ht node going down if position above had a child that implied an implicit node 
            # last condition to deal with case where we have moved down from the original graft point not on a speciation node,
            # in that case we prefer to look down 
            if(FALSE && currPos != finalPos && downHT[idNum[j]] && !nonAncRoot[idNum[j]] &&(currPos-numLeaves) > 0 && nodeType[currPos-numLeaves] == 'S'){ 
                #print('ht implicit check')
                implicitFlag = FALSE
                # get children of currPos
                childVec = tree$tree$edge[tree$tree$edge[,1] == currPos,2]
                # See if there is an inferred node at currPos
                #currPosID = speciesNamesToID[match(specIDall[currPos],speciesNamesToID[,3]),4]
                currPosID = specIDall[currPos]
                w = which(jpnames == currPosID)
                for(kc in 1:length(childVec)){
                    if(childVec[kc] >= rootpos && nodeType[childVec[kc]-numLeaves] != 'S') next
                    temp = speciesNamesToID[match(specIDall[childVec[kc]],speciesNamesToID[,3]),4]
                    while(!is.na(temp)){
                        if(temp %in% jpnames){
                            if(which(jpnames == temp) < w) implicitFlag = TRUE
                            break
                        }
                        temp = speciesNamesToID[match(temp,speciesNamesToID[,3]),4] 
                    }     
                }                               
                if(implicitFlag) {
                    implicitHT[idNum[j]] = TRUE
                    finalPos = currPos
                    if(finalPos != graftIDPos){
                        pathBetween = c(graftIDPos)
                        while(TRUE){
                            pathBetween = c(pathBetween,tree$tree$edge[tree$tree$edge[,2] == pathBetween[length(pathBetween)],1])
                            if(pathBetween[length(pathBetween)] == currPos) break
                        }
                        typeTemp = nodeType[pathBetween-numLeaves]
                        if(sum(typeTemp == 'D' ) > 0){
                            upDup[idNum[j]] = TRUE
                            dupUpList[[idNum[j]]] =  pathBetween[typeTemp == 'D']
                        }
                        if(sum(typeTemp == 'T' ) > 0){
                            upHT[idNum[j]] = TRUE
                            htUpList[[idNum[j]]] = pathBetween[typeTemp == 'T']
                        }
                    }
                
                    # Assigning other statisics 
                    testpos = 2
                    newAN[idNum[j]] = tree$tree$node.label[finalPos-numLeaves]
                    newGraft[idNum[j]] = finalPos
                    distBranches[idNum[j]] = tree2D[graftIDPos,finalPos]
                    distFloat[idNum[j]] = treeD[graftIDPos,finalPos]                 
                    graftNode = finalPos
                    #get pathway between the two graft points
                    currLoc = graftNode #note graftNode is position in node list of updated graft point
                    #pathway = pathBetween # inclusive of initial and final
                    graftCurr = graftIDPos
                     
                }
            }


            # If the new graft point ends up directly above a duplication node that was between it and the initial 
            # graft position it may be that there are implicit nodes closer to direction we came from : 
            #   1. For each other branch off the updated graft point
            #       A. Take all other branches from current node (recursively for duplication nodes) to get set of descendent speciations
            #       B. Check if there is an ancestor between them and current graft point to suggest there is implicit node in their direction
            #   2. Else move in opposite direction of original movement to the first speciation node, or if none, original graft point
            #       A. Unless there is an ht node in that direction
            # unless the other branches have an implicit node

            if(FALSE && lookup && finalPos %in% (parentVec+numLeaves)){ # if we are in a parent vector of original graft point
                #print('dup implicit check')
                m = match(finalPos,(parentVec+numLeaves))
                sChildren = c(tree$tree$edge[tree$tree$edge[,1]==finalPos,2])
                # Conditions : 1/2. we have to have moved up and there be a duplication below, 3. at least one speciation node as child 
                if( m > 1 && parentType[m-1] == 'D' && (sum(sChildren < rootpos) > 0 || sum(nodeType[sChildren-numLeaves] == 'S') > 0 ) ){ # if we moved up and there is a non speciation node right below 
                    # we will recurse through and get all child nodes not in parentVec that are speciation
                    w = which(sChildren %in% (parentVec+numLeaves))
                    if(length(w) > 1) print('error : to many child branches in parent vector ')
                    if(length(w) > 0) sChildren = sChildren[-w]
                    inferredState = FALSE
                    if(length(sChildren) < 1) print(paste('Error: no children that are speciation nodes but this shouldnt be possible',i,j))
                    for(kc in 1:length(sChildren)){
                        if(sChildren[kc] >= rootpos && nodeType[sChildren[kc]-numLeaves] != 'S') next
                        temp = speciesNamesToID[match(specIDall[sChildren[kc]],speciesNamesToID[,3]),4]
                        #templist = c(temp)
                        testpos = 3
                        finPosID = specIDall[finalPos]
                        w = which(jpnames == finPosID) #####################
                        while(!is.na(temp)){
                            temp = speciesNamesToID[match(temp,speciesNamesToID[,3]),4] 
                            #templist = c(templist,temp)
                            if(temp %in% jpnames){
                                if(which(jpnames == temp) < w) inferredState = TRUE
                                break
                            }
                        } 
                        if(inferredState) break   
                    }       


                    #for(ksc in 1:length(sChildren)){
                    #    while(TRUE){
                    #        taxtemp = speciesNamesToID[match(specIDall[sChildren],speciesNamesToID[,3]),4]
                    #        if(taxtemp %in% parentVecID || is.na(taxtemp)){
                    #            if(taxtemp == parentVecID[m]) interredState = TRUE
                    #            break
                    #        }
                    #    }
                    #    if(inferredState) break
                    #}
                    if(!inferredState){ # we will move if there is not an implicit node there
                        #dontMove = FALSE
                        inferredMove[idNum[j]] = TRUE
                        if(sum(parentType[1:(m-1)] != 'D') > 0 ) { # if we can move to a closer non dup do so 
                            w = which(parentType[1:(m-1)] == 'S')
                            w2 = which(parentType[1:(m-1)] == 'T')
                            # we will move, either to the first speciation or above the first ht
                            # this could put us on a duplication node in the second case
                            if(length(w) > 0 && length(w2) > 0 && max(w2) > max(w)) {
                                finalPos = parentVec[max(w2)+1] + numLeaves  #dontMove = TRUE
                            } else if(length(w) > 0) {
                                finalPos = parentVec[max(w)] + numLeaves 
                            } else if(length(w2) > 0 && length(w) == 0) {
                                finalPos = parentVec[max(w2)+1] + numLeaves 
                            }
                        }
                        else { # otherwise go back to graft point
                            finalPos = graftIDPos
                        }
                        if(finalPos >= rootpos && finalPos != graftIDPos){ # removed the dont move case

                            m2 = match(finalPos,parentVec+numLeaves)
                            w = m2 -1
                            testpos = 4
                            newAN[idNum[j]] = tree$tree$node.label[finalPos-numLeaves]
                            newGraft[idNum[j]] = finalPos
                            distBranches[idNum[j]] = tree2D[graftIDPos,finalPos]
                            distFloat[idNum[j]] = treeD[graftIDPos,finalPos]                 
                            graftNode = finalPos
                            #get pathway between the two graft points
                            currLoc = graftNode #note graftNode is position in node list of updated graft point
                            pathway = -1 # inclusive of initial and final
                            graftCurr = graftIDPos
                            if(w > 0){ 
                                parvec = parentVec[1:w][parentType[1:w] == 'D']
                                upDup[idNum[j]] = length(parvec) > 0
                                dupUpList[[idNum[j]]] = parvec
                            }
                            else {
                                upDup[idNum[j]] = FALSE
                                dupUpList[[idNum[j]]] = -1 # original graft point this should be nothing
                            }
                            
                        }
                    }


                } 


            }
             


            # If the new graft position is on the other side of a duplication node, and we did not pass through a ht node
            # move back down to first speciation node, or to original graft point if not
            if(FALSE && lookup && finalPos != graftIDPos && downDup[idNum[j]] && upDup[idNum[j]]){
                mr = mrcapairs(tree$tree$edge,finalPos,graftIDPos)
                if(nodeType[mr-numLeaves] == 'D' && !(downHT[idNum[j]]) && mr != rootpos && mr != graftIDPos){
                    m = match(mr-numLeaves,parentVec)
                    # bind to first speciation in parent list, or if there is an ht node, bind to first dup
                    ws = which(parentType[1:m] == 'S')
                    wh = which(parentType[1:m] == 'T')
                    #case1 = FALSE                    
                    if(length(ws) > 0 && length(wh) > 0 && max(wh) > max(ws)) {
                        finalPos = parentVec[max(wh)+1] + numLeaves
                        w = max(wh)
                    }
                    else if (length(ws) > 0) {
                        finalPos = parentVec[max(ws)] + numLeaves
                        w = max(ws)-1
                    }
                    else {
                        finalPos = graftIDPos
                        #case1 = TRUE
                        w = 0
                    }
                    switchBranch[idNum[j]] = TRUE
                    if(finalPos >= rootpos & finalPos != graftIDPos){
                        testpos = 5
                        newAN[idNum[j]] = tree$tree$node.label[finalPos-numLeaves]
                        newGraft[idNum[j]] = finalPos
                        distBranches[idNum[j]] = tree2D[graftIDPos,finalPos]
                        distFloat[idNum[j]] = treeD[graftIDPos,finalPos]                 
                        graftNode = finalPos
                        #get pathway between the two graft points
                        currLoc = graftNode #note graftNode is position in node list of updated graft point
                        pathway = -1 # inclusive of initial and final
                        graftCurr = graftIDPos
                        if(w > 0){ 
                            parvec = parentVec[1:w][parentType[1:w] == 'D']
                            upDup[idNum[j]] = length(parvec) > 0
                            dupUpList[[idNum[j]]] = parvec
                        }
                        else {
                            upDup[idNum[j]] = length(parvec) > 0
                            dupUpList[[idNum[j]]] = -1 # original graft point this should be nothing
                        }
                        
                        downDup[idNum[j]] = FALSE
                        dupDownList[[idNum[j]]] = -1
                    } 
                }
            }



            # If we ended up bound to a leaf due to movinig back to original graft point : 
            #   Set graft point right above leaf and add leaf to extra ortholog list

            if(finalPos < rootpos && finalPos == graftIDPos){
                #print('leaf fix')

                newAN[idNum[j]] = tree$tree$node.label[parentVec[1]]
                newGraft[idNum[j]] = parentVec[1]
                leafOrtho[idNum[j]] = graftIDPos
                distBranches[idNum[j]] = tree2D[graftIDPos,parentVec[1]+numLeaves]
                distFloat[idNum[j]] = treeD[graftIDPos,parentVec[1]+numLeaves]                 
                graftNode = parentVec[1]
                downDup[idNum[j]] = FALSE
                dupDownList[[idNum[j]]] = -1
                downHT[idNum[j]] = FALSE
                htDownList[[idNum[j]]] = -1
                upDup[idNum[j]] = FALSE
                dupUpList[[idNum[j]]] = -1
                upHT[idNum[j]] = FALSE
                htUpList[[idNum[j]]] = -1
                orthos[[idNum[j]]] = -7
                next

            }
            if(finalPos < rootpos && finalPos != graftIDPos){
                print(paste('Error : final position is a leaf but not graft point? ',i,j,finalPos,rootpos,graftIDPos))
                orthos[[idNum[j]]] = -6
                joasidjfasdf[[w3rfjo]]
                next
                #finalPos[[Inf]]
            }
            if(finalPos == graftIDPos){
                newGraft[idNum[j]] =  graftIDPos
                newAN[idNum[j]] = graftAN[j]
                pathways[[idNum[j]]] = -1
                currLoc = graftIDPos
                graftNode = graftIDPos
                distBranches[idNum[j]] = 0
                distFloat[idNum[j]] = 0
                orthos[[idNum[j]]] = '-0'
                next
            }     





            #### check if there are length 2 branches between initial and final graft point, these branches represent uncertain length
            if(distFloat[idNum[j]] >= 2){
        
                tree2.0 = tree
                tree2.0$tree$edge.length[treeDup$tree$edge.length == 2] = 100000
                treeD2.0 = dist.nodes(tree2.0$tree)
                dist2.0branches[idNum[j]] = floor(treeD2.0[newGraft[idNum[j]],graftIDPos] / 100000)

            }




            ##################################
            #######  Create Ortho Sets #######
            ##################################
     
            
            
            if(length(warnings()) > 0){
                print(warnings())
                 joaisdf[[Inf]]
            }
            if(newAN[idNum[j]] != an[idNum[j]]){ # should be true if it reaches this point
                o = orthosmrca[[match(unqpthr[i],names(orthosmrca))]]
                oposition = match(newAN[[idNum[j]]],tree$tree$node.label) # position in the list of node labels for graft
                o2 =  unpackorthos(o,oposition)  # gets list of positions in uniprot that will be orthos + xenos
                htg = htgroups[[match(unqpthr[i],names(htgroups))]] # get the htgroup
                htg2 = htgroups2[[match(unqpthr[i],names(htgroups))]] # get the htgroup branch version
                if(is.null(htg)) htg = rep(1,length(specIDall)) # dealing with case there are no ht groups
                if(is.null(htg2)) htg2 = rep(1,length(specIDall)) # dealing with case there are no ht groups
                hnum = htg[numLeaves + oposition]  
                hnum2 = htg2[numLeaves + oposition] 
                #htgOrthosNum = htg[hnum] What exactly was the point of this, think this is just outright incorrect
                hnumVec[i] = hnum
                htAmbig[[idNum[j]]] = c(uniprot[htg[1:numLeaves] == hnum & htg2[1:numLeaves] != hnum2]) 
                if(length( htAmbig[[idNum[j]]]) > 0) htAmbigtf[idNum[j]] = TRUE
                orthos[[idNum[j]]] = c(uniprot[o2 == 0 & htg2[1:numLeaves] == hnum2 & specleaves != sampSpec[j]]) # if we were bound to leaf after taking into account inferred nodes, add it back into set here
                if(leafOrtho[idNum[j]] != 0) orthos[[idNum[j]]] = c(orthos[[idNum[j]]],uniprot[leafOrtho[idNum[j]]])
                orthos[[idNum[j]]] = unique(orthos[[idNum[j]]]) # probably a faster way to do this
                xenos[[idNum[j]]] = uniprot[htg2[1:numLeaves] != hnum2]
                paras[[idNum[j]]] = uniprot[o2 == 1 & htg2[1:numLeaves] == hnum2] 
                if(FALSE && unkAmbig[[idNump[j]]]) unkAmbig[[idNum[j]]] = uniprot[htg2[1:numLeaves] != hnum2]
            } else{
                print('Error : At ortho creation but not different AN numbers?')
            }
            #children = get.descendants(tree,newGraft[idNum[j]])
            #children = children[children <= numLeaves]
            #childrenvec[idNum[j]] = str_c(uniprot[o2 == 0][children],collapse=';')
            #if(excludeChildren) orthos[[idNum[j]]] =  orthos[[idNum[j]]] = uniprot[o2 == 0][-children]




      

        }
        # TO DO keep below:? 
        if(FALSE & i == 100){ #for saving intermediate states
            print(warnings())
            print('saving intermediate state')
            uniprotList = rep('',i)
            for(counti in 1:length(orthos[pthr %in% unqpthr[1:i]])){
                uniprotList[counti] = str_c(orthos[pthr %in% unqpthr[1:i]][[counti]],collapse=';')
            }
            tempressamp = sampleName[pthr %in% unqpthr[1:i]]
            temprespthr = pthr[pthr %in% unqpthr[1:i]]
            tempresan = an[pthr %in% unqpthr[1:i]]
            tempnewan = newAN[pthr %in% unqpthr[1:i]]
            tempnewgraft = newGraft[pthr %in% unqpthr[1:i]]
            tempdistFloat = distFloat[pthr %in% unqpthr[1:i]]
            tempdistBranches = distBranches[pthr %in% unqpthr[1:i]]
            tempnote = note[pthr %in% unqpthr[1:i]]
            tempdownDup = downDup[pthr %in% unqpthr[1:i]]

            results = data.frame(tempressamp,temprespthr,tempresan,tempnewan,tempnewgraft,tempdistFloat,tempdistBranches,tempnote,tempdownDup,uniprotList)
            write.table(results,'./twoSpeciesTest11-1k.out',sep='\t')#str_c(saveloc,'loxaf-orthos-1500.out'),sep='\t')
        }
    }
    
#traceback()



#print(data.frame(sampleName,newAN,newGraft))   
uniprotList = rep('',length(pthr))
for(i in 1:length(orthos)){
    uniprotList[i] = str_c(orthos[[i]],collapse=';')
}
xenosList = rep('',length(pthr))
for(i in 1:length(orthos)){
    xenosList[i] = str_c(xenos[[i]],collapse=';')
}
 parasList = rep('',length(pthr))
for(i in 1:length(orthos)){
    parasList[i] = str_c(paras[[i]],collapse=';')
}


 
# could add in parasList and xenosList

results = data.frame(sampleName,inputSpecies,pthr,an,newAN,newGraft,distFloat,distBranches,jumpUpDist,downDup,upDup,downHT,upHT,inferredMove,implicitHT,tgRoot,leafGraft,nonAncRoot,matchBlockList,exceedBlockList,leafOrtho,switchBranch,note,uniprotList)
write.table(results,'./6speciesTest-2025-ver6-usingOldMoveDown-2.out',sep='\t')#str_c(saveloc,'loxaf-orthos-1500.out'),sep='\t')
# -2 was here before
print(paste('Output saved at : ',saveloc))

##################################################
################# Create Plot ####################
##################################################
if(makeplots){
    source('Code-Figures/generatePlot.r')
    for(i in 1:length(results[,1])) createFigure(r[i,])
}
        
  
# htgOrthosNumVec = rep(-1,length(pthr))
# hnumVec = rep(-1,length(pthr))