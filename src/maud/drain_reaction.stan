real drain_reaction(vector metabolite,
					real drain){
	fraction_of_drain = 1
	for (m in 1:size(metabolite)){
		fraction_of_drain *= metabolite[i] / (metabolite[i] + 1e-6)
	return rate * fraction_of_drain
}