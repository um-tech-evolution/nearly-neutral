#Rewriting  data/9_9_17Richness_plots.R in Julia to do plots using Plots package.
# 2/16/21
using Plots
using Measures
using CSV, Printf
using Colors, ColorTypes
#pwd()
function rich_hetero_IQV_plot( nmu::Float64, ptype::Symbol; savepng::Bool=false, unix::Bool=false,
      legend_position=:bottomleft, bg_color_legend=RGBA(1,1,1,1) )
  savedir = "3_22_21"   # subdirectory of data21 to save plots to.
  if !( nmu in [0.5,1.0,2.0,4.0] )
    error("nmu must be in [0.5,1.0,2.0,4.0]")
  end  
  if ptype in [:richness,:w_heteroz,:IQV]
    if !( nmu in [0.5,1.0,2.0,4.0] )
      error("nmu must be in [0.5,1.0,2.0,4.0]")
    end  
    (neut,del05,del10,mix05,mix1,mix2) = read_rich_hetero_IQV_data( unix=unix)
    println("data file read")
    selected_fields = [:N,:N_mu,:expected_w_heteroz,:w_heteroz,:expected_richness,:average_richness,:IQV]
    mu_sym = :N_mu
  elseif ptype == :num_fixed
    if !( nmu in [0.025,0.01,0.005,0.0025] )
      error("nmu must be in [0.025,0.01,0.005,0.0025]" ) 
    end  
    (neut,del05,del10,mix05,mix1,mix2) = read_num_mut_fixed( unix=unix)
    selected_fields = [:N,:mu,:num_extinct,:num_fixed,:fraction_fixed,:ave_fixed_time]
    mu_sym = :mu
  end
  
  neutstats = neut[neut[mu_sym].==nmu,selected_fields]
  del05stats=del05[del05[mu_sym].==nmu,selected_fields]
  del10stats=del10[del10[mu_sym].==nmu,selected_fields]
  mix05stats=mix05[mix05[mu_sym].==nmu,selected_fields]
  mix1stats=mix1[mix1[mu_sym].==nmu,selected_fields]
  mix2stats=mix2[mix2[mu_sym].==nmu,selected_fields]
  neut[!,:type]=fill(1,size(neut)[1])
  del05[!,:type]=fill(2,size(del05)[1])
  del10[!,:type]=fill(3,size(del10)[1])
  mix05[!,:type]=fill(4,size(mix05)[1])
  mix1[!,:type]=fill(5,size(mix1)[1])
  mix2[!,:type]=fill(6,size(mix2)[1])
  all_stats = DataFrames.vcat(neutstats, del05stats, del10stats,mix05stats,mix1stats,mix2stats)
  X = [25*2^i for i = 0:8]
  # Recreate xlabel and xticks since I can't get what I want through plots
  if ptype == :richness
    plottype = [:expected_richness,:average_richness]
  elseif ptype == :w_heteroz
    plottype = [:expected_w_heteroz,:w_heteroz]
  elseif ptype == :IQV
    plottype = [:IQV,:IQV]
  elseif ptype == :num_fixed
    plottype = [:num_fixed,:num_fixed]
  else
    error("plottype must be one of :richness, :w_heteroz, :IQV or :num_fixed")
  end
  plots = [
      neutstats[:,plottype[1]],
      neutstats[:,plottype[2]],
      del05stats[:,plottype[2]],
      del10stats[:,plottype[2]],
      mix05stats[:,plottype[2]],
      mix1stats[:,plottype[2]],
      mix2stats[:,plottype[2]]
  ]
  min_y = minimum([minimum(plots[i]) for i = 1:length(plots)])
  max_y = maximum([maximum(plots[i]) for i = 1:length(plots)])
  rng_y = max_y - min_y
  # Recreate xlabel and xticks since I can't get what I want through plots
  ano = vcat([(x,min_y-0.08*rng_y,text(@sprintf("%d",x),9)) for x in X],[(380,min_y-0.15*rng_y,text("Population size",10,))])
  labels = ["Expected" "Neutral" "Deleterious β=0.5" "Deleterious β=1.0" "Mixed β=0.005" "Mixed β=0.1" "Mixed β=0.2"]
  mcolors = [:brown :green :red :darkblue :coral4 :cyan :orange]
  mshapes = [:pentagon :rect :circle :utriangle :diamond :rtriangle :ltriangle]
  if ptype==:IQV
    plots = plots[2:end]
    ano = ano[2:end]
    labels = labels[:,2:end]
    mcolors = mcolors[:,2:end]
    mshapes = mshapes[:,2:end]
  end
  plot(X,plots,labels=labels,xscale=:log10,ann=ano,xticks=:none,bottom_margin=Measures.Length(:mm,12),color=mcolors,shape=mshapes)
  if ptype == :richness
    plot!(legend=:topleft,legendfontsize=9,bg_color_legend=bg_color_legend,title="Richness with θ=$(2*nmu)",ylabel="Richness",guidefontsize=12)
  elseif ptype == :w_heteroz
    plot!(legend=legend_position,legendfontsize=8,bg_color_legend=bg_color_legend,title="w heterozygosity with θ=$(2*nmu)",ylabel="Heterozygosity",guidefontsize=12)
  elseif ptype == :IQV
    plot!(legend=legend_position,bg_color_legend=bg_color_legend,legendfontsize=8,title="IQV with θ=$(2*nmu)",ylabel="IQV",guidefontsize=12)
  elseif ptype == :num_fixed
    plot!(legend=legend_position,bg_color_legend=bg_color_legend,legendfontsize=8,title="number mutations fixed with mu=$(nmu)",ylabel="Number fixed",guidefontsize=12)
  end
  if savepng
    cd(homedir())
    if unix
      cd("evotech/nearly-neutral")
    else
      cd("OneDrive/Research/NearlyNeutral")
    end
    cd("data21")
    savefig("$(savedir)/$(String(ptype))_theta=$(@sprintf("%3.1f",2*nmu)).png")
  end
end

# Read the csv files for the rich_hetero_IQV plots.

function read_rich_hetero_IQV_data(;unix::Bool=false)
  #datasubdir = "ngens200K"
  #datasubdir = "ngens400K"
  #datasubdir = "ngens800K"
  datasubdir = "9_9_17"
  datadir = "data17/"
  #datadir = "results21/Richness_figure3/"
  cd(homedir())
  if unix
    cd("evotech/nearly-neutral")
  else
    cd("OneDrive/Research/NearlyNeutral")
  end
  include("src/dataframe.jl")
  #cd("data17/9_9_17")
  cd("$(datadir)$(datasubdir)")
  println("read rich  pwd: ",pwd())
  neut = read_dataframe("nn_neut_Nmulist_fm45_N25_6400.csv")
  del05 = read_dataframe("nn_del_Nmulist_fm45_N25_6400_theta0.5.csv")
  del10 = read_dataframe("nn_del_Nmulist_fm45_N25_6400_theta1.csv")
  mix05 = read_dataframe("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.005.csv")
  mix1 = read_dataframe("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.01.csv")
  mix2 = read_dataframe("nn_mixed_Nmulist_fm45_prob0.01_N25_6400_theta0.5_0.02.csv")
  return (neut,del05,del10,mix05,mix1,mix2)
end

# Read the csv files for the num_mut_fixed plots.

function read_num_mut_fixed(;unix::Bool=false)
  cd(homedir())
  if unix
    cd("evotech/nearly-neutral")
  else
    cd("OneDrive/Research/NearlyNeutral")
  end
  include("src/dataframe.jl")
  cd("plots6_17/Figure4")
  println("pwd: ",pwd())
  neut = read_dataframe("in_neut_mulist_N25_6400.csv")
  del05 = read_dataframe("in_del_mulist_N25_6400_theta0.5.csv")
  del10 = read_dataframe("in_del_mulist_N25_6400_theta1.csv")
  mix05 = read_dataframe("in_mixed_mulist_N25_6400_theta0.5_0.005.csv")
  mix1 = read_dataframe("in_mixed_mulist_N25_6400_theta0.5_0.01.csv")
  mix2 = read_dataframe("in_mixed_mulist_N25_6400_theta0.5_0.02.csv")
  return (neut,del05,del10,mix05,mix1,mix2)
end

