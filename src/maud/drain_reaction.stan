real drain_reaction(vector metabolite,
					real drain){
  return drain * prod(metabolite ./ (metabolite + 1e-6));
}