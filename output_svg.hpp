#pragma once

/*SVG template
<svg
  id="mysvg"
  xmlns="http://www.w3.org/2000/svg"
  viewBox="0 0 800 600"
>
  <circle id="mycircle" cx="400" cy="300" r="50" />
</svg>
*/
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <utils.hpp>

class OutputSVG {
private:
	double _trans_base, _viewW, _viewH;
	std::string res;
	std::ofstream file;

	std::vector<double> trans(double x, double y) {
		std::vector<double> tmp = { 
			x - _trans_base, 
			_viewH - (y - _trans_base) 
		};
		return tmp;
	}

public:
	OutputSVG() {
		_trans_base = 0;
		_viewW = 1;
		_viewH = 1;
	}
	OutputSVG(double trans_base, double viewW, double viewH, std::string filename = "vor2d.svg") {
		_trans_base = trans_base;
		_viewW = viewW;
		_viewH = viewH;

		res = std::format(
			"<svg xmlns = \"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" width=\"500\" height=\"500\">\n",
			_viewW, _viewH);
	}

	void gen_line(double x1, double y1, double x2, double y2, std::string color = "#575757", double strokeWidth = 0.01) {
		//<line x1="0" y1="80" x2="100" y2="20" stroke="black" />
		auto p1 = trans(x1, y1),
			p2 = trans(x2, y2);
		std::string line = std::format(
			"<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"{}\"  />\n",
			p1[0], p1[1], p2[0], p2[1], color, strokeWidth);
		res += line;
	}
	void gen_point(double cx, double cy, std::string color = "red", double r = 0.001, double opacity = 1) {
		//<circle cx="50" cy="50" r="50" />
		auto p = trans(cx, cy);
		std::string circle = std::format(
			"<circle cx=\"{}\" cy=\"{}\" r=\"{}\" stroke=\"{}\" fill=\"none\" stroke-width=\"0.03\" fill-opacity=\"{}\"/>\n",
			p[0], p[1], r, color, opacity);
		res += circle;
	}
	void gen_any(std::string str) {
		res += str;
	}

	bool is_open() {
		return file.is_open();
	}
	void open_file(std::string filename) {
		file.open(UTILS::data_dir + filename);
		std::cout << UTILS::data_dir + filename << std::endl;
		res.clear();
		res = std::format(
			"<svg xmlns = \"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" width=\"500\" height=\"500\">\n",
			_viewW, _viewH);
	}
	void close_file() {
		res += "</svg>";
		file << res;
		file.close();
	}

};
