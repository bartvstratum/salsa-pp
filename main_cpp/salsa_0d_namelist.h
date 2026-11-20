/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Reads one group of a Fortran namelist file, for the subset used by the reference cases: one item per line,
 * values separated by commas, strings in quotes, logicals as .TRUE. or .FALSE., comments after '!'.
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>


namespace Salsa_0d
{
    // Item name in lower case, as Fortran is case insensitive, and its values as written in the file.
    using Namelist_group = std::map<std::string, std::vector<std::string>>;


    inline std::string trim(const std::string& text)
    {
        const size_t first = text.find_first_not_of(" \t\r");
        if (first == std::string::npos)
            return "";

        const size_t last = text.find_last_not_of(" \t\r");
        return text.substr(first, last - first + 1);
    }


    inline std::string to_lower(std::string text)
    {
        std::transform(text.begin(), text.end(), text.begin(), [](const unsigned char c) { return std::tolower(c); });
        return text;
    }


    inline std::string remove_comment(const std::string& line)
    {
        bool in_string = false;

        for (size_t i=0; i<line.size(); ++i)
        {
            if (line[i] == '\'' || line[i] == '"')
                in_string = !in_string;
            else if (line[i] == '!' && !in_string)
                return line.substr(0, i);
        }

        return line;
    }


    // Splits at the commas outside strings and removes the quotes around strings.
    inline std::vector<std::string> split_values(const std::string& text)
    {
        std::vector<std::string> values;
        std::string value;
        bool in_string = false;

        for (const char c : text)
        {
            if (c == '\'' || c == '"')
                in_string = !in_string;
            else if (c == ',' && !in_string)
            {
                values.push_back(trim(value));
                value.clear();
            }
            else
                value += c;
        }

        values.push_back(trim(value));
        return values;
    }


    inline Namelist_group read_namelist_group(const std::string& file_name, const std::string& group_name)
    {
        std::ifstream file(file_name);
        if (!file)
            throw std::runtime_error("Cannot open " + file_name);

        Namelist_group group;
        bool in_group = false;
        std::string line;

        while (std::getline(file, line))
        {
            line = trim(remove_comment(line));

            if (line.empty())
                continue;

            if (!in_group)
            {
                if (to_lower(line) == "&" + to_lower(group_name))
                    in_group = true;
                continue;
            }

            if (line == "/")
                return group;

            const size_t equals = line.find('=');
            if (equals == std::string::npos)
                throw std::runtime_error("Unsupported line in namelist group &" + group_name + ": " + line);

            const std::string name = to_lower(trim(line.substr(0, equals)));
            group[name] = split_values(line.substr(equals + 1));
        }

        throw std::runtime_error("Namelist group &" + group_name + " missing or not closed in " + file_name);
    }


    template<typename T>
    T convert(const std::string& value, const std::string& name)
    {
        try
        {
            if constexpr (std::is_same_v<T, bool>)
            {
                const std::string v = to_lower(value);
                if (v == ".true." || v == ".t." || v == "t")
                    return true;
                else if (v == ".false." || v == ".f." || v == "f")
                    return false;
                else
                    throw std::invalid_argument(value);
            }
            else if constexpr (std::is_same_v<T, std::string>)
                return value;
            else if constexpr (std::is_integral_v<T>)
                return std::stoi(value);
            else
            {
                // Fortran double precision exponents (1.d-3).
                std::string v = value;
                std::replace(v.begin(), v.end(), 'd', 'e');
                std::replace(v.begin(), v.end(), 'D', 'e');
                return static_cast<T>(std::stod(v));
            }
        }
        catch (const std::logic_error&)
        {
            throw std::runtime_error("Namelist item " + name + " has an invalid value: " + value);
        }
    }


    template<typename T>
    std::vector<T> get_list(const Namelist_group& group, const std::string& name)
    {
        const auto item = group.find(name);
        if (item == group.end())
            throw std::runtime_error("Namelist item " + name + " missing");

        std::vector<T> values;
        for (const std::string& value : item->second)
            values.push_back(convert<T>(value, name));

        return values;
    }


    template<typename T>
    T get_item(const Namelist_group& group, const std::string& name)
    {
        const std::vector<T> values = get_list<T>(group, name);
        if (values.size() != 1)
            throw std::runtime_error("Namelist item " + name + " must have one value");

        return values[0];
    }


    template<typename T>
    T get_item(const Namelist_group& group, const std::string& name, const T default_value)
    {
        if (group.find(name) == group.end())
            return default_value;

        return get_item<T>(group, name);
    }
}
