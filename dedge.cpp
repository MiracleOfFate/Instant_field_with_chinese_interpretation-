#include "dedge.h"

void build_dedge(const MatrixXu & F, const MatrixXf & V, VectorXu & V2E, VectorXu & E2E, VectorXb & boundary, VectorXb & nonManifold, const ProgressCallback & progress, bool quiet)
{
	// quiet = false ʱ
	if (!quiet) {
		// ������������ݽṹ
		cout << "Building a directed edge data structure .. ";
		cout.flush();
	}
	Timer<> timer;

	//progress = true �� quiet = false ʱ
	if (progress && !quiet)
		progress("Building directed edge data structure", 0.0f);

	// V2E ��ʼ����V2E��СΪ ������ * 1 ������
	V2E.resize(V.cols());
	V2E.setConstant(INVALID);

	uint32_t deg = F.rows();	// deg = 3
	//cout << "deg=" << deg << endl;
	std::vector<std::pair<uint32_t, uint32_t>> tmp(F.size());	// temp ��С = 3 * ����

	// ������
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)F.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t f = range.begin(); f != range.end(); ++f) {
			// һ�����������
			for (uint32_t i = 0; i < deg; ++i) {
				uint32_t idx_cur = F(i, f),			// ��ǰ��������
					idx_next = F((i + 1) % deg, f),	// ��һ��������
					edge_id = deg * f + i;			// ������

				// �������ݰ���Խ�綥��ο�
				if (idx_cur >= V.cols() || idx_next >= V.cols())
					throw std::runtime_error("Mesh data contains an out-of-bounds vertex reference!");
				
				// �Ի�
				if (idx_cur == idx_next)
					continue;

				tmp[edge_id] = std::make_pair(idx_next, INVALID);
				/*cout << "i=" << i << endl;
				cout << "edge_id=" << edge_id << endl;*/

				// Ϊ false ʱ
				if (!atomicCompareAndExchange(&V2E[idx_cur], edge_id, INVALID)) {
					uint32_t idx = V2E[idx_cur];
					//cout << "idx=" << idx << endl;
					while (!atomicCompareAndExchange(&tmp[idx].second, edge_id, INVALID)) {
						idx = tmp[idx].second;
						//cout << "tmp=" << tmp[idx].second << endl;
					}
				}
				//cout << "V2E=" << V2E << endl;
			}
		}
		if (!quiet)
			SHOW_PROGRESS_RANGE(range, F.cols(), "Building directed edge data structure (1/3)");
	}
	);
	/*cout << "V2E=" << V2E << endl;
	for (size_t i = 0; i < F.size(); i++)
	{
		cout << "tmp=" << tmp[i].first << "," << tmp[i].second << endl;
	}*/


	// nonManifold ��ʼ����nonManifold ��СΪ ������ * 1 ������
	nonManifold.resize(V.cols());
	nonManifold.setConstant(false);

	// E2E ��ʼ����E2E ��СΪ ������ * 3�� * 1 ������
	E2E.resize(F.cols() * deg);
	E2E.setConstant(INVALID);

	// ������
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)F.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t f = range.begin(); f != range.end(); ++f) {
			// һ�����������
			for (uint32_t i = 0; i < deg; ++i) {
				uint32_t idx_cur = F(i, f),			// ��ǰ��������
					idx_next = F((i + 1) % deg, f),	// ��һ��������
					edge_id_cur = deg * f + i;		// ������

				// �Ի�
				if (idx_cur == idx_next)
					continue;

				uint32_t it = V2E[idx_next], edge_id_opp = INVALID;
				/*cout << "it=" << it << endl;*/
				while (it != INVALID) {
					if (tmp[it].first == idx_cur) {
						if (edge_id_opp == INVALID) {
							edge_id_opp = it;
						}
						else {
							nonManifold[idx_cur] = true;
							nonManifold[idx_next] = true;
							edge_id_opp = INVALID;
							break;
						}
						/*cout << "edge_id_opp=" << edge_id_opp << endl;
						cout << "nonManifold=" << nonManifold << endl;*/
					}
					/*cout << "tmp[it].first=" << tmp[it].first << endl;
					cout << "idx_cur=" << idx_cur << endl;*/
					it = tmp[it].second;
				}

				if (edge_id_opp != INVALID && edge_id_cur < edge_id_opp) {
					E2E[edge_id_cur] = edge_id_opp;
					E2E[edge_id_opp] = edge_id_cur;
				}
			}
		}
		if (!quiet)
			SHOW_PROGRESS_RANGE(range, F.cols(), "Building directed edge data structure (2/3)");
	}
	);

	/*cout << "E2E=" << E2E << endl;
	cout << "nonManifold=" << nonManifold << endl;*/

	// �����ζ��������߽綥����������������
	std::atomic<uint32_t> nonManifoldCounter(0), boundaryCounter(0), isolatedCounter(0);

	boundary.resize(V.cols());
	boundary.setConstant(false);

	/* 
		Detect boundary regions of the mesh and adjust vertex->edge pointers 
		�������ı߽����򲢵��� ����->�� ָ��
	*/
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)V.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		// ��������
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			uint32_t edge = V2E[i];
			if (edge == INVALID) {
				isolatedCounter++;
				continue;
			}
			if (nonManifold[i]) {
				nonManifoldCounter++;
				V2E[i] = INVALID;
				continue;
			}

			/* 
				Walk backwards to the first boundary edge (if any)
				����ߵ���һ���߽�ߣ�����У�
			*/
			uint32_t start = edge, v2e = INVALID;
			do {
				v2e = std::min(v2e, edge);
				uint32_t prevEdge = E2E[dedge_prev(edge, deg)];
				if (prevEdge == INVALID) {
					/*
						Reached boundary -- update the vertex->edge link
						����߽硪�����¶���->�� ����
					*/
					v2e = edge;
					boundary[i] = true;
					boundaryCounter++;
					break;
				}
				edge = prevEdge;
			} while (edge != start);
			V2E[i] = v2e;
		}
		if (!quiet)
			SHOW_PROGRESS_RANGE(range, V.cols(), "Building directed edge data structure (3/3)");
	}
	);

	/*cout << "V2E=" << V2E << endl;
	cout << "E2E=" << E2E << endl;
	cout << "nonManifold=" << nonManifold << endl;
	cout << "boundary=" << boundary << endl;*/

	// quiet = false ʱ����ӡ��Ӧ�ķ����ζ��������߽綥����������������
	if (!quiet) {
		cout << "done. (";
		if (nonManifoldCounter)
			cout << nonManifoldCounter << " non-manifold vertices, ";
		if (boundaryCounter)
			cout << boundaryCounter << " boundary vertices, ";
		if (isolatedCounter)
			cout << isolatedCounter << " isolated vertices, ";
		cout << "took " << timeString(timer.value()) << ")" << endl;
	}
}