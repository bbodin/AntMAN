/*
 * testutils.h
 *
 *  Created on: May 11, 2021
 *      Author: toky
 */

#ifndef TESTS_CPP_TESTUTILS_H_
#define TESTS_CPP_TESTUTILS_H_

#include <boost/test/included/unit_test.hpp>
#include <AntMANLogger.h>
#include <utils.h>

std::string getString (AntMANLogger& logger) {
	 std::ostringstream outputstr;
	std::vector<std::string> names = {};
	int total = 0;

	if (logger.haslog("CI")) {total++;names.push_back("CI");}
	if (logger.haslog("CI")) {outputstr << "Size of CI vector = " << logger.getlog<cluster_indices_t>("CI").size() << std::endl ;}


	return outputstr.str();



}

#endif /* TESTS_CPP_TESTUTILS_H_ */
