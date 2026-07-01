#pragma once
#include <exception>
#include <string>

namespace stroid::exceptions {
    class StroidError : public std::exception {
    public:
        explicit StroidError(std::string message) : m_msg(std::move(message)) {}
        const char* what() const noexcept override { return m_msg.c_str(); }
    private:
        std::string m_msg;
    };

    class StroidContinuityError : public StroidError {
        using StroidError::StroidError;
    };

    class StroidMeshError : public StroidError {
        using StroidError::StroidError;
    };

    class StroidMissingReferenceMesh : public StroidMeshError {
        using StroidMeshError::StroidMeshError;
    };
}
