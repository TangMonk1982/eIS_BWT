#ifndef MY_TUPLE_H
#define MY_TUPLE_H

#include <iostream>
#include <limits>
#include <utility>  // for std::move

// 跨平台 packed 属性宏（与原代码行为一致）
#if defined(_MSC_VER)
    #define PACKED_BEGIN __pragma(pack(push, 1))
    #define PACKED_END   __pragma(pack(pop))
    #define PACKED_DECL
#elif defined(__GNUC__) || defined(__clang__)
    #define PACKED_BEGIN
    #define PACKED_END
    #define PACKED_DECL __attribute__((packed))
#else
    #define PACKED_BEGIN
    #define PACKED_END
    #define PACKED_DECL
#endif

//------------------------ Pair ------------------------//
PACKED_BEGIN
template<typename T1, typename T2>
struct PACKED_DECL Pair {
    T1 first;
    T2 second;

    // 构造函数
    Pair() : first(), second() {}                              // 默认构造（值初始化）
    Pair(const Pair&) = default;                               // 拷贝构造
    Pair(Pair&&) = default;                                    // 移动构造
    Pair(const T1& f, const T2& s) : first(f), second(s) {}    // 左值参数构造
    Pair(T1&& f, T2&& s) : first(std::move(f)), second(std::move(s)) {} // 右值参数构造

    // 赋值运算符（必须显式 default，否则因移动构造存在而被删除）
    Pair& operator=(const Pair&) = default;                    // 拷贝赋值
    Pair& operator=(Pair&&) = default;                          // 移动赋值

    // 最小值/最大值（返回非常量引用，与原代码兼容）
    static Pair& min_value() {
        static Pair min_val;
        return min_val;
    }
    static Pair& max_value() {
        static Pair max_val(std::numeric_limits<T1>::max(),
                            std::numeric_limits<T2>::max());
        return max_val;
    }

    // 相等比较
    bool operator==(const Pair& other) const {
        return first == other.first && second == other.second;
    }
    bool operator!=(const Pair& other) const {
        return !(*this == other);
    }

    // 输出运算符（修正为引用传递，不附加 endl）
    friend std::ostream& operator<<(std::ostream& os, const Pair& p) {
        os << "first: " << p.first << " second: " << p.second;
        return os;
    }
};
PACKED_END

//------------------------ Triple ------------------------//
PACKED_BEGIN
template<typename T1, typename T2, typename T3>
struct PACKED_DECL Triple {
    T1 first;
    T2 second;
    T3 third;

    Triple() : first(), second(), third() {}
    Triple(const Triple&) = default;
    Triple(Triple&&) = default;
    Triple(const T1& f, const T2& s, const T3& t) : first(f), second(s), third(t) {}
    Triple(T1&& f, T2&& s, T3&& t) : first(std::move(f)), second(std::move(s)), third(std::move(t)) {}

    // 赋值运算符
    Triple& operator=(const Triple&) = default;
    Triple& operator=(Triple&&) = default;

    static Triple& min_value() {
        static Triple min_val;
        return min_val;
    }
    static Triple& max_value() {
        static Triple max_val(std::numeric_limits<T1>::max(),
                              std::numeric_limits<T2>::max(),
                              std::numeric_limits<T3>::max());
        return max_val;
    }

    bool operator==(const Triple& other) const {
        return first == other.first && second == other.second && third == other.third;
    }
    bool operator!=(const Triple& other) const {
        return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const Triple& t) {
        os << "first: " << t.first << " second: " << t.second << " third: " << t.third;
        return os;
    }
};
PACKED_END

//------------------------ quadruple ------------------------//
PACKED_BEGIN
template<typename T1, typename T2, typename T3, typename T4>
struct PACKED_DECL quadruple {
    T1 first;
    T2 second;
    T3 third;
    T4 forth;  // 保留原拼写 forth

    quadruple() : first(), second(), third(), forth() {}
    quadruple(const quadruple&) = default;
    quadruple(quadruple&&) = default;
    quadruple(const T1& f, const T2& s, const T3& t, const T4& fo)
        : first(f), second(s), third(t), forth(fo) {}
    quadruple(T1&& f, T2&& s, T3&& t, T4&& fo)
        : first(std::move(f)), second(std::move(s)), third(std::move(t)), forth(std::move(fo)) {}

    // 赋值运算符
    quadruple& operator=(const quadruple&) = default;
    quadruple& operator=(quadruple&&) = default;

    static quadruple& min_value() {
        static quadruple min_val;
        return min_val;
    }
    static quadruple& max_value() {
        static quadruple max_val(std::numeric_limits<T1>::max(),
                                 std::numeric_limits<T2>::max(),
                                 std::numeric_limits<T3>::max(),
                                 std::numeric_limits<T4>::max());
        return max_val;
    }

    bool operator==(const quadruple& other) const {
        return first == other.first && second == other.second &&
               third == other.third && forth == other.forth;
    }
    bool operator!=(const quadruple& other) const {
        return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const quadruple& q) {
        os << "first: " << q.first << " second: " << q.second
           << " third: " << q.third << " forth: " << q.forth;
        return os;
    }
};
PACKED_END

//------------------------ quintuple ------------------------//
PACKED_BEGIN
template<typename T1, typename T2, typename T3, typename T4, typename T5>
struct PACKED_DECL quintuple {
    T1 first;
    T2 second;
    T3 third;
    T4 forth;
    T5 fifth;

    quintuple() : first(), second(), third(), forth(), fifth() {}
    quintuple(const quintuple&) = default;
    quintuple(quintuple&&) = default;
    quintuple(const T1& f, const T2& s, const T3& t, const T4& fo, const T5& fi)
        : first(f), second(s), third(t), forth(fo), fifth(fi) {}
    quintuple(T1&& f, T2&& s, T3&& t, T4&& fo, T5&& fi)
        : first(std::move(f)), second(std::move(s)), third(std::move(t)),
          forth(std::move(fo)), fifth(std::move(fi)) {}

    // 赋值运算符
    quintuple& operator=(const quintuple&) = default;
    quintuple& operator=(quintuple&&) = default;

    static quintuple& min_value() {
        static quintuple min_val;
        return min_val;
    }
    static quintuple& max_value() {
        static quintuple max_val(std::numeric_limits<T1>::max(),
                                 std::numeric_limits<T2>::max(),
                                 std::numeric_limits<T3>::max(),
                                 std::numeric_limits<T4>::max(),
                                 std::numeric_limits<T5>::max());
        return max_val;
    }

    bool operator==(const quintuple& other) const {
        return first == other.first && second == other.second &&
               third == other.third && forth == other.forth && fifth == other.fifth;
    }
    bool operator!=(const quintuple& other) const {
        return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const quintuple& q) {
        os << "first: " << q.first << " second: " << q.second
           << " third: " << q.third << " forth: " << q.forth
           << " fifth: " << q.fifth;
        return os;
    }
};
PACKED_END

#endif // MY_TUPLE_H
