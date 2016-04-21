#ifndef FILEIO_H

#define FILEIO_H

#include "std.h"
#include "data.h"
#include "io.h"
#include "functions.h"

class fileio{ 
	public:
        static void read_map (string filename, unordered_map<string,int> *m, int idx = 0, int vidx = 0);
        static void read_map (string filename, unordered_map<string,double> *m, int idx = 0, int vidx = 0);
		static void read_map (string file, unordered_map<string,string> *m, int idx = 0 , int vidx = -1);
        static void read_vector (string filename, vector<string> *m, int idx = -1);


        template<typename T>
            static std::istream & binary_read(std::istream& stream, T& value){
                return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
            } 

        template<typename T>
            static std::istream & binary_read(std::istream& stream, vector<T>& v){
                int size;
                binary_read (stream, size);
                v.resize (size);
                for (int i = 0 ; i < v.size(); i++) {
                    binary_read (stream, v[i]);
                }

                return stream;
            } 

        template<typename T>
            static std::ostream& binary_write(std::ostream& stream, const T& value){
                return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
            }

        template<typename T>
            static std::ostream& binary_write(std::ostream& stream, const vector<T>& v){
                int size = v.size();
                binary_write (stream, size);
                for (int i = 0 ; i < v.size(); i++)
                    binary_write (stream, v[i]);
                return stream;

            }


};
#endif
