#pragma once

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <stack>
#include <set>
#include <cstdlib>
#include <cassert>
#include <Eigen/Core>
#include <cgal_types_2d.hpp>
#include <output_svg.hpp>


namespace TEST_2D {
	class PowerDiagram {
	private:
		Regular_triangulation rt;
		Eigen::MatrixXd Vertices, Normals;
		Eigen::VectorXd Weights;
		int PSize;

		double q_offset = -0.2;  //fixed epsilon

		Iso_rectangle_2 bbox = Iso_rectangle_2(Point(-4, -4), Point(4, 4));


	private:
		void gen_weights(double w_diff, double w_base = 4);
		void generate_verts(double d);
		void compute_regular_triangulation();
		bool crop_and_extract_segment(CGAL::Object& o, Segment_2& res_seg);
		void get_voronoi(
			Eigen::MatrixXd& vertices,  //n*2
			Eigen::MatrixXi& segments_idx  //m*2
		);


	public:
		void core_main(double d);



	};

}