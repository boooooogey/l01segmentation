#include "range.hpp"
#include <iterator>

Range::Range(int capacity){
    this->capacity = 2 * capacity;
    data = new double[this->capacity];
    length = 0;
}

Range::Range(const Range& other){
    if(this != &other){
        length = other.length;
        capacity = other.capacity;
        data = new double[capacity];
        std::copy(other.data, other.data + length, data);
    }
}

int Range::len() const{
    return length;
}

void Range::add(const double& first, const double& last){
    if(length == capacity){
        capacity = 2 * capacity;

        double* tmpdata = data;
        data = new double[capacity];
        std::copy(tmpdata, tmpdata+length, data);
        delete[] tmpdata;
    }
    data[2*length] = first;
    data[2*length + 1] = last;
    length++;
}

void Range::index(const int& i, double& first, double& last) const{
    if(length-1 < i){
        return;
    }
    else{
        first = data[2*i];
        last = data[2*i + 1];
    }
}

Range::~Range(){
    delete[] data;
}

RangeList::RangeList(){
    length = 0;
    data = nullptr;
}

RangeList::RangeList(const int& length){
    this->length = length;
    data = new Range[this->length];
}

RangeList::RangeList(const RangeList& other){
    if(this != &other){
        length = other.length;
        delete[] data;
        data = new Range[length];
        std::copy(other.data, other.data + length, data);
    }
}

void RangeList::resize(const int& length){
    if(data != nullptr){
        delete[] data;
    }
    this->length = length;
    data = new Range[this->length];
}

int RangeList::len() const{
    return length;
}

int RangeList::len(const int& i) const{
    return data[i].len();
}

void RangeList::add(const int& i, const double& first, const double& last){
    if(0 <= i && i < length)
        data[i].add(first, last);
}

void RangeList::index(const int& i, const int& j, double& first, double& last) const{
    if(length-1 < i || i < 0){
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
    delete[] data;
}
