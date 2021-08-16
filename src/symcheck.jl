# Read in a config file, write out the structures one by one, check symmetry
using DelimitedFiles

"""
Parse an array of lines (from an MTP-style config file) into a structure. Also return the unused lines. This routine is fragile in the sense that it assumes spacing between consecutive structures remains constant.
"""
function parsecfg!(lines)
popfirst!(lines) # Skip BEGIN_CFG
popfirst!(lines) # Skip Size
Nats = parse(Int64,popfirst!(lines))
popfirst!(lines) # Skip Supercell
LV = zeros(3,3)
LV[1,:] = parse.(Float64,split(popfirst!(lines)))
LV[2,:] = parse.(Float64,split(popfirst!(lines)))
LV[3,:] = parse.(Float64,split(popfirst!(lines)))
popfirst!(lines) # Skip AtomData:  id type...
aTyp = zeros(Int64,Nats)
pos = zeros(Nats,3)
for i ∈ 1:Nats
    aTyp[i] = parse.(Int64,split(lines[1])[2])
    pos[i,:] = parse.(Float64,split(lines[1])[3:5])
    popfirst!(lines)
end
while !contains(popfirst!(lines),"END_CFG") end # skip ahead to end of configuration
popfirst!(lines) # Skip blank line
return LV, aTyp, pos
end

function writeSymcheckFile(lattvecs,atypes,positions)
open("/Users/hart/codes/symlib/src/symcheck.str","w") do f
    writedlm(f,lattvecs)
    writedlm(f,length(atypes))
    for i ∈ 1:length(atypes)
        writedlm(f,Any[atypes[i] positions[i,:]'])
    end
end
end

function getSPcount(out)
nSPops = []
while length(out) > 1
    lattvecs, atypes, positions = parsecfg!(out)
    writeSymcheckFile(lattvecs,atypes,positions)
    symx=`/Users/hart/codes/symlib/src/sym.x /Users/hart/codes/symlib/src/symcheck.str`
    push!(nSPops,parse(Int64,split(read(symx,String))[1]))
end
return nSPops
end

path=pwd()*"/preselected.cfg" 
out = open(path) do file 
    readlines(file) end 
splist=getSPcount(out)

# pull the conf_id out of cfg file and parse the string to get the size of the pointgroup
cfg = open(path) do file 
    readlines(file) end 
cflist=[parse(Int64,split(i,"_")[3][3:end]) for i ∈ filter(p->contains(p,"conf_id"),cfg)]
# This next line handles cfg files with the newer format (due to updated makeStr.py output)
append!(cflist,[parse(Int64,split(i)[3]) for i ∈ filter(p->contains(p,"sym_size"),cfg)])
# Get the starting and ending line numbers of each structure in the cfg file
startidx=findall(x->x=="BEGIN_CFG",cfg)
endidx  =findall(x->x=="END_CFG",cfg)


# Loop over each structures and separate into two files, those that preserved the symmetry
# and those where symmetry was broken
safefile=pwd()*"/safesympreselected.cfg"
brokenfile=pwd()*"/brokenpreselect.cfg"
rm(safefile,force=true); rm(brokenfile,force=true)
for (i,str) in enumerate(splist)
    println("Checking str#: ",i)
    if cflist[i]==splist[i] # symmetry has been preserved, print this to the training file.
        println("Found a valid structure")
        targetfile = safefile
    else
        println("Found a broken structure")
        targetfile =  brokenfile
    end
    open(targetfile, "a") do io
        for line in cfg[startidx[i]:endidx[i]] 
            println(io, line) 
        end
    end
end

