#include "range.hpp"
#include <iterator>

Range::Range(int capacity){
    data.reserve(2 * capacity);
    data.resize(0);
}

Range::Range(const Range& other){
    if(this != &other){
        data = other.data;
    }
}

int Range::len() const{
    return data.size()/2;
}

void Range::add(const double& first, const double& last){
    data.push_back(first);
    data.push_back(last);
}

void Range::index(const int& i, double& first, double& last) const{
    if(data.size()-1 < i){
        return;
    }
    else{
        first = data[2*i];
        last = data[2*i + 1];
    }
}

Range::~Range(){
}

RangeList::RangeList(){
    data.resize(0);
}

RangeList::RangeList(const int& length){
    data.reserve(length);
    data.resize(length);
}

RangeList::RangeList(const RangeList& other){
    if(this != &other){
        data = other.data;
    }
}

void RangeList::resize(const int& length){
    data.resize(length);
}

int RangeList::len() const{
    return data.size(); 
}

int RangeList::len(const int& i) const{
    return data[i].len();
}

void RangeList::add(const int& i, const double& first, const double& last){
    if(0 <= i && i < data.size())
        data[i].add(first, last);
}

void RangeList::index(const int& i, const int& j, double& first, double& last) const{
    if(data.size()-1 < i || i < 0){
        return;
    }
    else{
        data[i].index(j, first, last);
    }
}

Range& RangeList::operator[](const int& i){
    return data[i];
}

RangeList::~RangeList(){
}
