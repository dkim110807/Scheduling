//
// Created by dkim110807 on 2024-06-11.
//

#pragma once

#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <bitset>

std::string to_string(const std::string &s) {
    return '"' + s + '"';
}

std::string to_string(const char *s) {
    return to_string((std::string) s);
}

std::string to_string(bool b) {
    return (b ? "true" : "false");
}

std::string to_string(int x) {
    return std::to_string(x);
}

std::string to_string(int64_t x) {
    return std::to_string(x);
}

std::string to_string(size_t x) {
    return std::to_string(x);
}

std::string to_string(char x) {
    std::string s;
    s += x;
    return s;
}

std::string to_string(const std::vector<bool> &v) {
    bool first = true;
    std::string res = "{";
    for (auto &&i: v) {
        if (!first) res += ", ";
        first = false;
        res += to_string(i);
    }
    res += "}";
    return res;
}

template<typename T, size_t N>
std::string to_string(const std::array<T, N> &v) {
    std::string res = "(";
    for (size_t i = 0; i < N; i++) {
        if (i > 0) res += ", ";
        res += to_string(v[i]);
    }
    res += ')';
    return res;
}

template<typename A>
std::string to_string(A v) {
    bool first = true;
    std::string res = "{";
    for (const auto &x: v) {
        if (!first) {
            res += ", ";
        }
        first = false;
        res += to_string(x);
    }
    res += "}";
    return res;
}

void debug_out() { std::cerr << std::endl; }

template<typename Head, typename... Tail>
void debug_out(Head H, Tail... T) {
    std::cerr << " " << to_string(H);
    debug_out(T...);
}

#ifdef LOCAL
#define debug(...) std::cerr << "Line " <<  __LINE__ << " => "<<  "[" << #__VA_ARGS__ << "] :", debug_out(__VA_ARGS__)
#else
#define debug(...) 0x110807
#endif