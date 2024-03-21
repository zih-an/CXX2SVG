#include <power_diagram.hpp>
using namespace TEST_2D;


void PowerDiagram::gen_weights(double w_diff, double w_base) {
	Weights.conservativeResize(Vertices.rows());
	for (int i = 0; i < Vertices.rows(); i++) {
		if (i < PSize)  Weights(i) = w_base + w_diff;
		else Weights(i) = w_base;
	}
}


void PowerDiagram::generate_verts(double d) {
	std::ifstream ifs(UTILS::data_dir + "points.txt");
	Vertices.conservativeResize(10, 2);
	Normals.conservativeResize(10, 2);

	std::string s;
	int cnt = 0;
	while (std::getline(ifs, s)) {
		std::istringstream str(s);
		double v1, v2, n1, n2;
		str >> v1 >> v2 >> n1 >> n2;
		Vertices.row(cnt) << v1, v2;
		Normals.row(cnt) << n1, n2;

		cnt++;
	}
	Vertices.conservativeResize(cnt, 2);
	Normals.conservativeResize(cnt, 2);
	Normals.rowwise().normalize();

/*----------------------------------------------------------------------------*/
	PSize = Vertices.rows();
	Vertices.conservativeResize(PSize * 2, Vertices.cols());
	Vertices.block(PSize, 0, PSize, Vertices.cols()) =
		Vertices.block(0, 0, PSize, Vertices.cols()) + q_offset * Normals;

	std::cout << Vertices << std::endl;
	std::cout << Normals << std::endl;

	double w_diff = 2 * q_offset * d - q_offset * q_offset;
	gen_weights(w_diff);
}


void PowerDiagram::compute_regular_triangulation() {
	rt.clear();
	std::vector<std::pair<Weighted_point, unsigned>> w_points;
	for (int i = 0; i < Vertices.rows(); i++) {
		w_points.push_back(std::make_pair(
			Weighted_point(
				Point(Vertices(i, 0), Vertices(i, 1)),
				Weights(i)
			),
			i
		));
	}

	rt.insert(w_points.begin(), w_points.end());
	assert(rt.is_valid());
}


bool PowerDiagram::crop_and_extract_segment(CGAL::Object& o, Segment_2& res_seg) {
	CGAL::Object obj;

	if (CGAL::object_cast<Segment_2>(&o)) {
		//std::cout << "Segment_2" << std::endl;
		if (CGAL::object_cast<Segment_2>(o).source() == CGAL::object_cast<Segment_2>(o).target())  //dual may be a point, but as seg with same end-points
			return false;
		obj = CGAL::intersection(CGAL::object_cast<Segment_2>(o), bbox);
	}
	else if (CGAL::object_cast<Ray_2>(&o)) {
		//std::cout << "Ray_2" << std::endl;
		obj = CGAL::intersection(CGAL::object_cast<Ray_2>(o), bbox);
	}
	else if (CGAL::object_cast<Line_2>(&o)) {
		//std::cout << "Line_2" << std::endl;
		obj = CGAL::intersection(CGAL::object_cast<Line_2>(o), bbox);
	}
	else {
		std::cout << "[wrong] oops!" << std::endl;
	}
	const Segment_2* s = CGAL::object_cast<Segment_2>(&obj);
	if (s != nullptr) {
		res_seg = *s;
		return true;
	}
	return false;
}

void PowerDiagram::get_voronoi(
	Eigen::MatrixXd& vertices,  //n*2
	Eigen::MatrixXi& segments_idx  //m*2
) {
	std::vector<double> edge_verts;  int cnt_verts = 0;
	std::vector<int> segs_idx;
	//int ns = 0, nr = 0;
	for (Edge_iterator eit = rt.edges_begin();
		eit != rt.edges_end(); ++eit
		) {
		Vertex_handle
			v1 = eit->first->vertex(Regular_triangulation::ccw(eit->second)),
			v2 = eit->first->vertex(Regular_triangulation::cw(eit->second));
		if ((v1->info() < PSize && v2->info() < PSize) || (v1->info() >= PSize && v2->info() >= PSize))
			continue;
		              

		CGAL::Object o = rt.dual(eit);
		Segment_2 seg;
		if (!crop_and_extract_segment(o, seg)) continue;


		edge_verts.push_back(CGAL::to_double(seg.source().x()));
		edge_verts.push_back(CGAL::to_double(seg.source().y()));
		//edge_verts.push_back(0);//z=0

		edge_verts.push_back(CGAL::to_double(seg.target().x()));
		edge_verts.push_back(CGAL::to_double(seg.target().y()));
		//edge_verts.push_back(0);//z=0

		segs_idx.push_back(cnt_verts++);
		segs_idx.push_back(cnt_verts++);

	}

	Eigen::MatrixXd tmp_verts = Eigen::Map<Eigen::MatrixXd>(edge_verts.data(), 2, cnt_verts);
	vertices = tmp_verts.transpose();
	Eigen::MatrixXi tmp_segs = Eigen::Map<Eigen::MatrixXi>(segs_idx.data(), 2, cnt_verts / 2);
	segments_idx = tmp_segs.transpose();

}


void PowerDiagram::core_main(double d) {
	std::cout << "generate verts" << std::endl;
	generate_verts(d);
	std::cout << "compute_regular_triangulation" << std::endl;
	compute_regular_triangulation();


//output
	Eigen::MatrixXd vor_verts;
	Eigen::MatrixXi vor_segs;
	std::cout << "get_voronoi" << std::endl;
	get_voronoi(vor_verts, vor_segs);

/*-----------------------------------------------------------------------------------*/
	OutputSVG svg;
	double xmin = CGAL::to_double(bbox.xmin()),
		xmax = CGAL::to_double(bbox.xmax()),
		ymin = CGAL::to_double(bbox.ymin()),
		ymax = CGAL::to_double(bbox.ymax());
	svg = OutputSVG(
		std::min(xmin, ymin),
		xmax - xmin,
		ymax - ymin
	);


	svg.open_file("result.svg");
	svg.gen_any("<g>\n");
	for (int i = 0; i < Vertices.rows(); i++) {
		Eigen::VectorXd p = Vertices.row(i);
		if (i<PSize) {
			svg.gen_point(p(0), p(1), "#ffb399", 0.1, 1);
			svg.gen_point(p(0), p(1), "#ffb399", std::sqrt(Weights(i)), 0.1);//die
		}
		else {
			svg.gen_point(p(0), p(1), "#c7ebff", 0.1, 1);
			svg.gen_point(p(0), p(1), "#c7ebff", std::sqrt(Weights(i)), 0.1);//die
		}

	}
	svg.gen_any("</g>\n");

	svg.gen_any("<g>\n");
	for (int i = 0; i < vor_segs.rows(); i++) {
		Eigen::VectorXd p1 = vor_verts.row(vor_segs(i, 0)),
			p2 = vor_verts.row(vor_segs(i, 1));
		svg.gen_line(p1(0), p1(1), p2(0), p2(1), "black", 0.01);
	}
	svg.gen_any("</g>\n");
	svg.close_file();

}
