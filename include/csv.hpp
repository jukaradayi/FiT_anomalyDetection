// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree
#include "history_graph.hpp"

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class CSVRow
{
    public:
        std::string operator[](std::size_t index) const
        {
            return std::string(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(' ', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   

class CSVIterator
{   
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) {++(*this); }
        CSVIterator()                   :m_str(NULL) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};

class CSVRange
{
    std::istream&   stream;
    public:
        CSVRange(std::istream& str)
            : stream(str)
        {}
        CSVIterator begin() const {return CSVIterator{stream};}
        CSVIterator end()   const {return CSVIterator{};}
};
class InputParser{

    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }

        const std::pair<std::vector<StreamGraphs::Count>, StreamGraphs::Count> getSizes() const{
            std::vector<StreamGraphs::Count> hist_sizes;
            std::vector<std::string>::const_iterator itrH =  std::find(this->tokens.begin(), this->tokens.end(), "-h");
            std::vector<std::string>::const_iterator itrG =  std::find(itrH, this->tokens.end(), "-g");
            std::vector<std::string>::const_iterator next_cmd = std::find(itrG+1, this->tokens.end(), "-k");

            StreamGraphs::Count G_index = 0;
            if (itrH != this->tokens.end()){
                for (++itrH; itrH != itrG; ++itrH) {
                    ++G_index;
                    const std::string& h_size = *itrH;
                    std::cout << h_size << "\n";

                    hist_sizes.push_back(std::stoi(h_size));
                }
            }

            if (itrG != this->tokens.end()){
                for (++itrG; itrG != next_cmd; ++itrG) {
    
                    const std::string& g_size = *itrG;
                    std::cout << g_size << "\n";

                    hist_sizes.push_back(std::stoi(g_size));
                }
            }
            std::pair<std::vector<StreamGraphs::Count>, StreamGraphs::Count> Sizes = std::make_pair(hist_sizes, G_index);
            return Sizes;
         
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

