#pragma once

#include <algorithm>
#include <boost/utility/string_ref.hpp>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace graphtools
{
// Borrowed text; conversion to an owning string is explicit at API boundaries
// through the standard string conversion, never during mapped catalog loading.
class TextView : public boost::string_ref
{
public:
    TextView() = default;
    TextView(const char* value)
        : boost::string_ref(value)
    {
    }
    TextView(const std::string& value)
        : boost::string_ref(value)
    {
    }
    TextView(const char* value, size_t length)
        : boost::string_ref(value, length)
    {
    }
    operator std::string() const { return std::string(data(), size()); }
    TextView substr(size_t start, size_t count = boost::string_ref::npos) const
    {
        const auto value = boost::string_ref::substr(start, count);
        return TextView(value.data(), value.size());
    }
};
inline bool operator==(TextView a, TextView b) { return boost::string_ref(a) == boost::string_ref(b); }
inline bool operator!=(TextView a, TextView b) { return !(a == b); }
inline bool operator==(TextView a, const char* b) { return a == TextView(b); }
inline bool operator==(const char* a, TextView b) { return b == a; }
inline bool operator!=(TextView a, const char* b) { return !(a == b); }
inline bool operator!=(const char* a, TextView b) { return !(b == a); }
inline std::string operator+(const std::string& a, TextView b) { return a + std::string(b); }
inline std::string operator+(const char* a, TextView b) { return std::string(a) + std::string(b); }
inline std::string& operator+=(std::string& a, TextView b) { return a.append(b.data(), b.size()); }
inline bool operator==(TextView a, const std::string& b) { return boost::string_ref(a) == boost::string_ref(b); }
inline bool operator==(const std::string& a, TextView b) { return b == a; }
inline bool operator!=(TextView a, const std::string& b) { return !(a == b); }
inline bool operator!=(const std::string& a, TextView b) { return !(a == b); }
inline std::ostream& operator<<(std::ostream& out, TextView text) { return out.write(text.data(), text.size()); }

struct Slice
{
    uint64_t offset, count;
};

inline uint64_t hashText(TextView text)
{
    uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : text)
    {
        hash ^= c;
        hash *= 1099511628211ULL;
    }
    return hash;
}

struct MappedBytes
{
    MappedBytes(const char* bytes = nullptr, size_t length = 0, std::shared_ptr<const void> lifetime = {})
        : data(bytes)
        , size(length)
        , owner(std::move(lifetime))
    {
    }
    const char* data = nullptr;
    size_t size = 0;
    std::shared_ptr<const void> owner;
    template <typename T> const T* array(Slice slice) const
    {
        if (slice.offset > size || slice.count > (size - slice.offset) / sizeof(T) || slice.offset % alignof(T))
            throw std::runtime_error("Invalid mapped array range or alignment");
        return reinterpret_cast<const T*>(data + slice.offset);
    }
    TextView text(Slice slice) const { return TextView(array<char>(slice), slice.count); }
};

class BlobBuilder
{
public:
    template <typename T> Slice append(const T* values, size_t count)
    {
        const size_t padding = (alignof(T) - bytes.size() % alignof(T)) % alignof(T);
        bytes.append(padding, '\0');
        Slice slice { bytes.size(), count };
        if (count)
            bytes.append(reinterpret_cast<const char*>(values), sizeof(T) * count);
        return slice;
    }
    template <typename T> Slice append(const T& value) { return append(&value, 1); }
    template <typename T> Slice append(const std::vector<T>& values) { return append(values.data(), values.size()); }
    Slice text(TextView value) { return append(value.data(), value.size()); }
    template <typename T> void replace(uint64_t offset, const T& value)
    {
        if (offset > bytes.size() || sizeof(T) > bytes.size() - offset)
            throw std::out_of_range("Blob patch");
        std::memcpy(&bytes[offset], &value, sizeof(T));
    }
    std::string bytes;
};

// Copyable range descriptors and iterators contain only process-local handles;
// the elements stay in the mapped file. Iterators own their descriptor so a
// range returned by value cannot leave a dangling iterator.
template <typename T> class CatalogRange
{
public:
    using Getter = T (*)(const void*, size_t, const void*, size_t);
    CatalogRange() = default;
    CatalogRange(
        const void* context, size_t size, Getter get, const void* aux = nullptr, size_t extent = 0,
        std::shared_ptr<const void> owner = {})
        : context_(context)
        , aux_(aux)
        , size_(size)
        , extent_(extent)
        , get_(get)
        , owner_(std::move(owner))
    {
    }
    template <typename U>
    static auto retain(U& value, const std::shared_ptr<const void>& owner, int)
        -> decltype(value.retainOwner(owner), void())
    {
        value.retainOwner(owner);
    }
    template <typename U> static void retain(U&, const std::shared_ptr<const void>&, long) { }
    class iterator
    {
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = T;
        iterator(
            const void* context, size_t index, Getter get, const void* aux, size_t extent,
            std::shared_ptr<const void> owner)
            : context_(context)
            , aux_(aux)
            , index_(index)
            , extent_(extent)
            , get_(get)
            , owner_(std::move(owner))
        {
        }
        T operator*() const
        {
            T value = get_(context_, index_, aux_, extent_);
            retain(value, owner_, 0);
            return value;
        }
        iterator& operator++()
        {
            ++index_;
            return *this;
        }
        iterator operator++(int)
        {
            auto old = *this;
            ++*this;
            return old;
        }
        iterator& operator--()
        {
            --index_;
            return *this;
        }
        iterator operator--(int)
        {
            auto old = *this;
            --*this;
            return old;
        }
        bool operator==(const iterator& other) const { return context_ == other.context_ && index_ == other.index_; }
        bool operator!=(const iterator& other) const { return !(*this == other); }

    private:
        const void* context_;
        const void* aux_;
        size_t index_;
        size_t extent_;
        Getter get_;
        std::shared_ptr<const void> owner_;
    };
    iterator begin() const { return iterator(context_, 0, get_, aux_, extent_, owner_); }
    iterator end() const { return iterator(context_, size_, get_, aux_, extent_, owner_); }
    std::reverse_iterator<iterator> rbegin() const { return std::reverse_iterator<iterator>(end()); }
    std::reverse_iterator<iterator> rend() const { return std::reverse_iterator<iterator>(begin()); }
    size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }
    T operator[](size_t index) const
    {
        T value = get_(context_, index, aux_, extent_);
        retain(value, owner_, 0);
        return value;
    }
    T front() const { return (*this)[0]; }
    T back() const { return (*this)[size_ - 1]; }
    operator std::vector<T>() const { return std::vector<T>(begin(), end()); }
    size_t count(const T& value) const { return std::count(begin(), end(), value); }
    iterator find(const T& value) const { return std::find(begin(), end(), value); }
    bool operator==(const CatalogRange& other) const
    {
        return size() == other.size() && std::equal(begin(), end(), other.begin());
    }
    bool operator!=(const CatalogRange& other) const { return !(*this == other); }
    static CatalogRange fromVector(const std::vector<T>& values, std::shared_ptr<const void> owner = {})
    {
        return CatalogRange(
            &values, values.size(), [](const void* p, size_t i, const void*, size_t)
            { return (*static_cast<const std::vector<T>*>(p))[i]; }, nullptr, 0, std::move(owner));
    }
    static CatalogRange fromArray(const T* values, size_t size, std::shared_ptr<const void> owner = {})
    {
        return CatalogRange(
            values, size, [](const void* p, size_t i, const void*, size_t) { return static_cast<const T*>(p)[i]; },
            nullptr, 0, std::move(owner));
    }
    template <typename Container>
    static CatalogRange fromContainer(const Container& values, std::shared_ptr<const void> owner = {})
    {
        return CatalogRange(
            &values, values.size(),
            [](const void* p, size_t i, const void*, size_t) -> T
            {
                auto it = static_cast<const Container*>(p)->begin();
                std::advance(it, i);
                return *it;
            },
            nullptr, 0, std::move(owner));
    }

private:
    const void* context_ = nullptr;
    const void* aux_ = nullptr;
    size_t size_ = 0;
    size_t extent_ = 0;
    Getter get_ = nullptr;
    std::shared_ptr<const void> owner_;
};

template <typename T, typename U, typename Alloc>
bool operator==(const std::vector<U, Alloc>& values, const CatalogRange<T>& view)
{
    return values.size() == view.size() && std::equal(values.begin(), values.end(), view.begin());
}
template <typename T, typename U, typename Compare, typename Alloc>
bool operator==(const std::set<U, Compare, Alloc>& values, const CatalogRange<T>& view)
{
    return values.size() == view.size() && std::equal(values.begin(), values.end(), view.begin());
}
template <typename T, typename U, typename Hash, typename Equal, typename Alloc>
bool operator==(const std::unordered_set<U, Hash, Equal, Alloc>& values, const CatalogRange<T>& view)
{
    if (values.size() != view.size())
        return false;
    for (const auto& value : view)
        if (!values.count(value))
            return false;
    return true;
}
}
