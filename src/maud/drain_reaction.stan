real drain_reaction(vector metabolite,
					real drain){
	real fraction_of_drain = 1;
	for (m in 1:size(metabolite)){
		fraction_of_drain *= metabolite[m] / (metabolite[m] + 1e-6);
  }
	return drain * fraction_of_drain;
}