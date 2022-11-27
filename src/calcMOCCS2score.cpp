#include <Rcpp.h>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <regex>
#include <thread>
#include <chrono>

using namespace Rcpp;

std::int8_t getKmerDistance(
    std::string fasta,
    std::unordered_map<std::string, std::vector<std::int32_t>> *uom,
    std::int32_t klen,
    std::int32_t *fastaLen,
    std::int32_t *ind,
    std::int32_t *seqLenMax
);

std::vector<std::double_t> calcAUCandMOCCS2score(
    std::vector<std::int32_t> posVec,
    std::int32_t seqLenMax,
    std::int32_t klen
);

//List incProgress(std::double_t oneStep, std::string message) {
//    Function f("incProgress");
//    return f(oneStep, Named("message") = message);
//}
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<std::double_t>> cppCalcMOCCS2score(
    std::string fileName,
    std::int32_t klen,
    bool ignoreLowerCase
) {
    std::unordered_map<std::string, std::vector<std::int32_t>> uom;
    std::unordered_map<std::string, std::vector<std::double_t>> retm;
    std::string fasta;
    std::int32_t count=0;
    std::int32_t lineNum=1;
    std::string strForCount;
    std::ifstream ifs(fileName);
    std::ifstream forCount(fileName);
    std::string fastaHead;
    std::int32_t fastaLen=0;
    std::int32_t ind=0;
    std::int32_t seqLenMax = 0;
    List hoge;

    if (ifs.fail()) {
        return retm;
    }
    while (getline(forCount, strForCount)) {
        ++count;
    }
    count = floor(count / 20);
    // std::double_t onestep = 1 / count;
    while (getline(ifs, fasta)) {
        if (lineNum >= count) {
            //incProgress(0.05, "Processing fasta file...");
            lineNum = 1;
            // std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
        fastaHead = fasta[0];
        if (fastaHead != ">") {
            if (ignoreLowerCase) {
                fasta = std::regex_replace(fasta, std::regex("[atgcn]"), "N");
            } else {
                std::transform(fasta.begin(), fasta.end(), fasta.begin(), toupper);
            }
            getKmerDistance(fasta, &uom, klen, &fastaLen, &ind, &seqLenMax);
        }
        lineNum++;
    }
  //  incProgress(0, "Caluculating MOCCS2score...");
    for (auto it = uom.begin(); it != uom.end(); ++it) {
        retm[it->first] = calcAUCandMOCCS2score(it->second, seqLenMax, klen);
    }
    return retm;
}

std::int8_t getKmerDistance(
    std::string fasta,
    std::unordered_map<std::string, std::vector<std::int32_t>> *uom,
    std::int32_t klen,
    std::int32_t *fastaLen,
    std::int32_t *ind,
    std::int32_t *seqLenMax
) {
    *fastaLen = fasta.length();
    *seqLenMax = std::max({*fastaLen, *seqLenMax});
    std::int32_t n = *fastaLen -klen + 1;
    std::double_t getCenter =(klen - *fastaLen) *0.5;
    std::string query;
    for(*ind = 0;*ind < n; ++*ind) {
        query = (std::string)fasta.substr(*ind, klen);
        if (query.find("N") == std::string::npos) {
            (*uom)[query].push_back(floor(abs(*ind + getCenter)));
        }
    }

    return 0;
}


std::vector<std::double_t> calcAUCandMOCCS2score(
    std::vector<std::int32_t> posVec,
    std::int32_t seqLenMax,
    std::int32_t klen
) {
    auto posvecbegin = posVec.begin();
    auto posvecend = posVec.end();
    std::int32_t vmax = floor(abs((seqLenMax - klen) * 0.5));
    std::vector<std::int32_t> countVec(vmax + 1);
    std::int32_t ind;
    for(ind=0; ind<vmax + 1; ++ind) {
        countVec[ind] = std::count(posvecbegin, posvecend, ind) * (vmax - ind);
    }
    std::double_t retVal = std::accumulate(countVec.begin(), countVec.end(), 0);
    std::int32_t posSize = posVec.size();
    retVal = (retVal / posSize) - (vmax * 0.5);
    std::double_t window = floor(seqLenMax - 1)/2 + 1 - floor(klen / 2);
    std::double_t moccs2Score = retVal / window * sqrt(12) * sqrt(posSize);
    std::vector<std::double_t> retvec = {
        retVal,
        static_cast<std::double_t>(posSize),
        moccs2Score,
        window
    };
    return retvec;
}
