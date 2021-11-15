#pragma once
#include <vector>

class Range{
    private:
        std::vector<double> data;
    public:
        Range(int capacity = 5);
        Range(const Range& other);
        int len() const;
        void add(const double& first, const double& last);
        void index(const int& i, double& first, double& last) const;
        ~Range();
};

class RangeList{
    private:
        std::vector<Range> data;
        int length;
    public:
        RangeList();
        RangeList(const int& length);
        RangeList(const RangeList& other);
        void resize(const int& length);
        int len() const;
        int len(const int& i) const;
        void add(const int& i, const double& first, const double& last);
        void index(const int& i, const int& j, double& first, double& last) const;
        Range& operator[](const int& i);
        ~RangeList();
};
