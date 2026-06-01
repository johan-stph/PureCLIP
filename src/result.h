// ======================================================================
// PureCLIP — lightweight Result<T> for error propagation
// ======================================================================
// Replaces boolean return + std::cerr side-effects with structured errors.
// Usage:
//   Result<void> r = doSomething();
//   if (!r) { std::cerr << r.error() << std::endl; return false; }
// ======================================================================

#ifndef PURECLIP_RESULT_H_
#define PURECLIP_RESULT_H_

#include <string>
#include <utility>
#include <variant>

template <typename T = void>
class Result {
    std::variant<T, std::string> storage_;
    bool ok_;

public:
    // Success constructors
    Result(T value) : storage_(std::move(value)), ok_(true) {}
    
    // Error constructor
    Result(std::string error) : storage_(std::move(error)), ok_(false) {}

    explicit operator bool() const noexcept { return ok_; }
    bool has_value() const noexcept { return ok_; }

    T& value() {
        return std::get<T>(storage_);
    }
    const T& value() const {
        return std::get<T>(storage_);
    }

    const std::string& error() const {
        return std::get<std::string>(storage_);
    }

    T value_or(T&& default_val) const {
        return ok_ ? std::get<T>(storage_) : std::forward<T>(default_val);
    }
};

// Specialization for void (no value on success, just error string on failure)
template <>
class Result<void> {
    std::string error_;

public:
    Result() = default;  // success
    Result(std::string error) : error_(std::move(error)) {}

    explicit operator bool() const noexcept { return error_.empty(); }
    bool has_value() const noexcept { return error_.empty(); }

    const std::string& error() const { return error_; }
};

#endif // PURECLIP_RESULT_H_
