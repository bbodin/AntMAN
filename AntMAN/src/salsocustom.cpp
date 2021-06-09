#include "salsocustom.h"

double negative_infinity = -std::numeric_limits<double>::infinity();

salso_result_t salsoCpp(const arma::mat &eam, ind_t maxClusts,
						double Const_Binder, ind_t batchSize, ind_t nScans,
						unsigned int maxThreads, unsigned int timeLimitPerThread)
{

	arma::mat p = eam - Const_Binder; // we rarely use eam directly
	p.diag().zeros();				  // items do not contribute to the Binder score by being in
	// the same cluster with themselves
	ind_t N = p.n_cols; // number of items
	salso_result_t result(N);
	if (maxClusts == 0) {
		maxClusts = N;
	} else {
		maxClusts = std::min(maxClusts, N); // never need more than N clusters
	}

	int numThreads;
#pragma omp parallel
	{
#pragma omp single
		{
			VERBOSE_DEBUG ( "Begin clustering " );
			numThreads = 1;
			result.numThreads = numThreads;
			if (numThreads == 1) {
				VERBOSE_DEBUG ( "using 1 thread.\n");
			} else {
				VERBOSE_DEBUG ( "using " << numThreads << " threads.\n");
			}
			VERBOSE_DEBUG ( "Number of permutations to search: " << numThreads * batchSize << '\n');
		}
		auto timeStart = std::chrono::high_resolution_clock::now();

		salso_result_t partialResult(N);	 // partial result from each thread
		labelVec_t cl(N, arma::fill::zeros); // cluster label vectors
		// iterate over random orderings of of [1, ..., N]
		while (true) {
			/* BEGIN ITERATION */
			double currIterBinderScore = 0;			  // stores the score for this iteration
			arma::uvec itemOrder(randperm(N));		  // Generate the random item ordering for this iteration. We will assign cluster labels to items in this order.
			arma::mat pOrd = p(itemOrder, itemOrder); // Permute the co-clustering probability matrix to match the item ordering for this iteration.
			/* FIRST SEQUENTIAL ALLOCATION
				Sequentially allocate labels to all N items as follows:
				1. Assign the 0th item to the 0th cluster.
				2. For k = 1, ..., N-1, given a clustering of the first k-1 items:
				a. Create all possible clusterings of the first k items by
				varying the kth item label only. b. The kth item label can be any of
				the labels already present, or a new label. b. Find the best among
				these clusterings.
			*/
			std::fill(cl.begin(), cl.end(), 0);						 // cl stores our current best labelling
			std::vector<std::vector<ind_t>> labelIndices(maxClusts); // labelIndices[t] stores the indices of the items with label t
			labelIndices[0].push_back(0);							 // item 0 gets label 0 to begin with
			ind_t currNumClusts = 1, tryNumClusts;					 // currently only using 1 cluster
			/* FIRST SEQUENTIAL ALLOCATION */
			for (ind_t k = 1; k < N; ++k) {
				/* ITERATE OVER ITEMS */
				tryNumClusts = std::min(currNumClusts + 1, maxClusts);
				double bestLabelDelta = negative_infinity, tmpLabelDelta; // change in Binder score
				ind_t bestLabel = 0;
				for (ind_t t = 0; t < tryNumClusts; ++t) {
					/* ITERATE OVER CANDIDATE LABELS */
					tmpLabelDelta = arma::accu(pOrd.unsafe_col(k).elem(arma::uvec(labelIndices[t]))); // change in Binder score if we use label t for item k
					if (tmpLabelDelta > bestLabelDelta) { // want to maximise the change
						bestLabelDelta = tmpLabelDelta;
						bestLabel = t;
					}
					/* END CURRENT LABEL */
				}
				labelIndices[bestLabel].push_back(k); // item k has label bestLabel
				cl[k] = bestLabel;					  // item k has label bestLabel
				if (bestLabel == tryNumClusts - 1 && currNumClusts < maxClusts) {
					currNumClusts++; // item was assigned a label not currently in our set of labels
				}
				currIterBinderScore += bestLabelDelta;
				/* END CURRENT ITEM */
			}
			/* END FIRST SEQUENTIAL ALLOCATION */

			/* SWEETENING SCANS
				1. Perform step 2 nScans times.
				2. For k = 0, ..., N-1:
				a. Consider all possible clusterings of the N items obtained
				by varying the kth item label only. b. The kth item label can be any
				of the labels already present, or a new label. c. Find the best
				among these clusterings.
				3. At the end of a scan if there is no change in clustering over the
				previous scan, stop scanning.
			*/
			ind_t currMaxLabel = currNumClusts - 1, tryMaxLabel; // stores the current maximum label and the label to try during the scan
			double thisScanDeltaBinder = 0;						 // change in Binder score from this scan
			for (ind_t currScan = 0; currScan < nScans; ++currScan) {
				/* BEGIN SCAN */
				for (ind_t k = 0; k < N; ++k) {
					/* ITERATE OVER ITEMS */
					tryMaxLabel = std::min(currMaxLabel + 1, maxClusts - 1);							   // ensures that we have at most maxClusts clusters. Leads to some redundant iterations if there are holes in the labelling, but does not make a huge difference to speed and is easier than bookkeeping the holes.
					double bestLabelDelta = 0, tmpLabelDelta, currLabelDelta;							   // worst we can do is nothing
					currLabelDelta = arma::accu(pOrd.unsafe_col(k).elem(arma::uvec(labelIndices[cl[k]]))); // contribution to current Binder score due to the kth item
					ind_t bestLabel = cl[k];
					for (ind_t t = 0; t <= tryMaxLabel; ++t) {
						/* ITERATE OVER CANDIDATE LABELS */
						tmpLabelDelta = -currLabelDelta + arma::accu(pOrd.unsafe_col(k).elem(arma::uvec(labelIndices[t]))); // change in Binder score if we change to label t for item k
						if (tmpLabelDelta > bestLabelDelta) { // want to maximise the change
							bestLabelDelta = tmpLabelDelta;
							bestLabel = t;
						}
						/* END CURRENT LABEL */
					}
					if (bestLabel == cl[k]) {
						continue; // no change in label
					}
					ind_t oldLabel = cl[k];
					labelIndices[oldLabel].erase(std::find(labelIndices[oldLabel].begin(), labelIndices[oldLabel].end(), k)); // remove item k from the list of items having its current label
					labelIndices[bestLabel].push_back(k);																	  // add item k to the list of items having label bestLabel
					cl[k] = bestLabel;																						  // item k has label bestLabel
					if (bestLabel > currMaxLabel)
						++currMaxLabel;
					thisScanDeltaBinder += bestLabelDelta;
					/* END CURRENT ITEM */
				}
				if (thisScanDeltaBinder == 0) {
					break; // no change in Binder score from the scan
				}
				currIterBinderScore += thisScanDeltaBinder;
				/* END SCAN */
			}

			if (currIterBinderScore > partialResult.binderLoss)
			{ // if the current iteration yielded a better clustering
#pragma omp simd
				for (ind_t k = 0; k < N; ++k) {
					partialResult.labels[itemOrder[k]] = cl[k]; // undo the permutation on the labels
				}
				partialResult.binderLoss = currIterBinderScore;
				/*partialResult.numClusts = cl.max ();*/
			}

			/* LOOP EXIT CONDITIONS AND BOOKKEEPING */
			++partialResult.nIters;
			auto timeNow = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(timeNow - timeStart).count();

			/* CHECK ITERATION COUNT */
			if (partialResult.nIters >= batchSize && batchSize > 0)
			{
				partialResult.wallClockTime = duration;
				if (timeLimitPerThread > 0 && duration >= timeLimitPerThread) {
					partialResult.timeLimitReached = true;
				}
				break;
			}

			/* CHECK TIMEOUT */
			if (timeLimitPerThread > 0 && duration >= timeLimitPerThread)
			{
				partialResult.wallClockTime = duration;
				partialResult.timeLimitReached = true;
				break;
			}

			/* END ITERATION */
		}

#pragma omp critical
		{
			result.nIters += partialResult.nIters;
			result.wallClockTime += partialResult.wallClockTime;
			result.timeLimitReached |= partialResult.timeLimitReached;
			if (partialResult.binderLoss > result.binderLoss)
			{
				result.labels = std::move(partialResult.labels);
				//result.numClusts = partialResult.numClusts;
				result.binderLoss = partialResult.binderLoss;
			}
		}
	}

	/* CANONICALISE LABELS STARTING AT 1 */
	labelVec_t sortedLabels(N, arma::fill::zeros), labelPerm(result.labels.max() + 1, arma::fill::zeros);
	for (ind_t i = 0, c = 0; i < N; ++i)
	{
		if (labelPerm[result.labels[i]] == 0) {
			labelPerm[result.labels[i]] = ++c;
		}
	}
#pragma omp simd
	for (ind_t i = 0; i < N; ++i){
		sortedLabels[i] = labelPerm[result.labels[i]];
	}
	result.labels = std::move(sortedLabels);
	result.numClusts = result.labels.max();

	/* ADJUST THE BINDER LOSS TO CORRECTLY REFLECT THE SCORE */
	result.binderLoss = -result.binderLoss + (1 - Const_Binder) * arma::accu(arma::trimatu(eam, 1));

	/* OUTPUT IF NOT IN R */
#ifndef HAS_RCPP
	VERBOSE_DEBUG ( "Cluster labels:");
	VERBOSE_DEBUG ( result.labels);
	VERBOSE_DEBUG ( "Finished clustering, found " << result.numClusts << " clusters.");
	VERBOSE_DEBUG ( "Normalised binder loss: " << result.binderLoss );
	VERBOSE_DEBUG ( "Number of permutations scanned: " << result.nIters);
	VERBOSE_DEBUG ( "Time limit reached: " << result.timeLimitReached);
#endif
	return result;
}

double computeBinderLossCpp(const arma::mat &eam,
							const arma::ivec &partitionLabels,
							double Const_Binder)
{
	arma::mat p = eam - Const_Binder;
	ind_t N = partitionLabels.n_elem;
	double Binder_f = 0;
	for (ind_t j = 0; j < N - 1; ++j)
	{
#pragma omp simd
		for (ind_t k = j + 1; k < N; ++k)
		{
			if (partitionLabels[j] == partitionLabels[k]) {
				Binder_f += p(j, k);
			}
		}
	}
	return -Binder_f + (1 - Const_Binder) * arma::accu(arma::trimatu(eam, 1));
}

std::vector<ind_t> randperm(ind_t N)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::vector<ind_t> ans(N);
	std::iota(ans.begin(), ans.end(), 0); // sequentially fill values
	//std::shuffle (ans.begin (), ans.end (), mt);
	return ans;
}
