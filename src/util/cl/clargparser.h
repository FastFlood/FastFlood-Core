#ifndef FULLFLOOD_CLARGPARSER_H
#define FULLFLOOD_CLARGPARSER_H

#include <vector>
#include "../util/stringutil.h"
#include <string>
#include "../util/log.h"
#include "../common/fferror.h"
class ArgParser
{

public:

    constexpr static std::string ARGP_TYPE_ANY ="ANY";
    constexpr static std::string ARGP_TYPE_DOUBLE ="DOUBLE";
    constexpr static std::string ARGP_TYPE_STRING ="STRING";
    constexpr static std::string ARGP_TYPE_BOOL ="BOOL";

private:
    std::vector<std::string> m_Names;
    std::vector<std::vector<std::string>> m_Values;
    std::vector<std::vector<std::vector<std::string>>> m_ValuesAdditive;
    std::vector<std::string> m_Description;
    std::vector<int> m_Lengths;
    std::vector<bool> m_Set;
    std::vector<bool> m_IsAdditive;
    std::vector<std::vector<std::string>> m_Types;
    std::vector<std::vector<std::string>> m_SubDescriptions;
    std::vector<bool> m_Required;

public:

    inline void AddOption(std::string Name, std::string description,  int length, std::vector<std::string> types,std::vector<std::string> subdescription, bool is_required, bool is_additive = false)
    {
        m_Names.push_back(Name);
        m_Description.push_back(description);
        m_Lengths.push_back(length);
        m_Types.push_back(types);
        m_SubDescriptions.push_back(subdescription);
        m_Required.push_back(is_required);
        m_IsAdditive.push_back(is_additive);
        m_Set.push_back(false);
        std::vector<std::string> valueslist ={};
        for(int i = 0; i < types.size(); i++)
        {
            valueslist.push_back("");
        }

        m_Values.push_back(valueslist);
        m_ValuesAdditive.push_back({});
    }
    inline void PrintOptionsAll()
    {

        PrintOptions(true);

    }
    inline void PrintOptions(bool all = false)
    {
        FF::coutI <<"Input options for FastFlood " << std::endl;
        for(int i = 0; i < m_Names.size(); i++)
        {
            FF::coutI << m_Names[i] + (m_IsAdditive[i]? " - (multiple supported) - " :" -  ") + m_Description[i] << std::endl;
            if(all)
            {
                FF::coutI << "Number of arguments: " << m_Values[i].size() << std::endl;
                for(int k = 0; k < m_Values[i].size(); k++)
                {
                    std::string req = " ";
                    if(m_Required[i])
                    {
                        req = " [REQUIRED]";
                    }
                    if(k > m_Types[i].size()-1 || k > m_SubDescriptions[i].size() - 1)
                    {
                        FF::coutI<<"CONFIG ERROR, PLEASE CONTACT FASTFLOOD" <<std::endl;
                    }else
                    {
                        FF::coutI << m_Types[i][k] << " - " <<  m_SubDescriptions[i][k] << req << std::endl;
                    }
                }
                FF::coutI << std::endl;
            }
        }
    }

    inline int Parse(std::string option, std::vector<std::string> args, bool force_overwrite = false)
    {
            //get argument name
            std::string str_name(option);
            FF_INFO("Check arg "+ str_name);
            if(str_name[0] == '-')
            {
                bool found = false;
                for(int j = 0; j < m_Names.size(); j++)
                {
                    if(str_name == ("-" + m_Names[j]))
                    {
                        if(!(!m_Set[j] || force_overwrite))
                        {
                            return -1;
                        }

                        FF_INFO("Set arg " + m_Names[j]);
                        found = true;
                        //set the flag as set
                        m_Set[j] = true;
                        //and get all the sub values
                        int a = 0;
                        if(m_IsAdditive[j])
                        {
                            m_ValuesAdditive[j].push_back({});
                            a = m_ValuesAdditive[j].size()-1;
                            for(int k = 0; k < m_Values[j].size(); k++)
                            {
                                m_ValuesAdditive[j][a].push_back("");
                            }
                        }

                            for(int k = 0; k < m_Values[j].size(); k++)
                            {

                                if(k < args.size())
                                {
                                    //FF::cout <<"GetVal For " << str_name << " " << k << " " << i+1+k << " " << argc << std::endl;
                                    //FF::cout << std::string(argv[i+1 + k]) <<std::endl;
                                    m_Values[j][k] = args[k];
                                    if(m_IsAdditive[j])
                                    {
                                        m_ValuesAdditive[j][a][k] = args[k];
                                    }

                                    FF_INFO("Arg value " + m_Values[j][k]);
                                    if(m_Types[j][k] == "DOUBLE")
                                    {
                                        if(!is_number(m_Values[j][k]))
                                        {
                                            ff_error = "Variable is not a number: " +m_Values[j][k] +" for " + m_Names[j];
                                            return -1;
                                        }

                                    }else if(m_Types[j][k] == "INT")
                                    {
                                        if(!is_number(m_Values[j][k]))
                                        {
                                            ff_error = "Variable is not a number: " + m_Values[j][k] + " for " + m_Names[j];
                                            return -1;
                                        }
                                    }else if(m_Types[j][k] == "STRING")
                                    {

                                    }
                                }else
                                {
                                    FF::coutE << "Expected " << m_Values[j].size() << " arguments for option " << str_name << " but did not find enough" << std::endl;
                                    return -1;
                                }
                            }



                        break;

                    }
                }
                if(!found)
                {
                    FF::coutE << "Parser could not parse unkown argument " << str_name << std::endl;
                }
            }

        return 0;


    }

    inline int Parse(int argc, const char ** argv)
    {
        for(int i = 0; i < argc; i++)
        {
            //get argument name
            std::string str_name(argv[i]);
            FF_INFO("Check arg "+ str_name);
            if(str_name[0] == '-')
            {
                bool found = false;
                for(int j = 0; j < m_Names.size(); j++)
                {
                    if(str_name == ("-" + m_Names[j]))
                    {
                        FF_INFO("Set arg " + m_Names[j]);
                        found = true;
                        //set the flag as set
                        m_Set[j] = true;
                        //and get all the sub values
                        int a = 0;

                        if(m_IsAdditive[j])
                        {
                            m_ValuesAdditive[j].push_back({});
                            a = m_ValuesAdditive[j].size()-1;
                            for(int k = 0; k < m_Values[j].size(); k++)
                            {
                                m_ValuesAdditive[j][a].push_back("");
                            }
                        }

                        if(argc > i + m_Values[j].size())
                        {
                            for(int k = 0; k < m_Values[j].size(); k++)
                            {
                                //FF::cout <<"GetVal For " << str_name << " " << k << " " << i+1+k << " " << argc << std::endl;
                                //FF::cout << std::string(argv[i+1 + k]) <<std::endl;
                                m_Values[j][k] = trimc(std::string(argv[i+1 + k]));
                                if(m_IsAdditive[j])
                                {
                                    m_ValuesAdditive[j][a][k] = trimc(std::string(argv[i+1 + k]));
                                }
                                if(m_Names[j] != "key" && m_Names[j] != "token")
                                {
                                    FF_INFO("Arg value " + m_Values[j][k]);
                                }
                                if(m_Types[j][k] == "DOUBLE")
                                {
                                    if(!is_number(m_Values[j][k]))
                                    {
                                        ff_error = "Variable is not a number: " +m_Values[j][k] +" for " + m_Names[j];
                                        return -1;
                                    }

                                }else if(m_Types[j][k] == "INT")
                                {
                                    if(!is_number(m_Values[j][k]))
                                    {
                                        ff_error = "Variable is not a number: " + m_Values[j][k] + " for " + m_Names[j];
                                        return -1;
                                    }
                                }else if(m_Types[j][k] == "STRING")
                                {

                                }
                            }

                            i = i +m_Values[j].size();

                        }else
                        {
                            FF::coutE << "Expected " << m_Values[j].size() << " arguments for option " << str_name << " but did not find enough" << std::endl;
                            return -1;
                        }

                        break;

                    }
                }
                if(!found)
                {
                    FF::coutE << "Parser could not parse unkown argument " << str_name << std::endl;
                }
            }
        }
        return 0;
    }
    inline void RemoveArg(std::string name)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                m_Set[j] =false;
            }
        }
    }

    inline bool IsNumber(std::string name, int i)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                return is_number(m_Values[j][i]);
            }
        }
        return false;
    }
    inline bool IsArgSet(std::string name)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                return m_Set[j];
            }
        }
        return false;
    }

    //specific for additive variables
    inline int GetArgSetCount(std::string name)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                if(m_IsAdditive[j])
                {
                    return m_ValuesAdditive[j].size();
                }else
                {
                    return m_Set[j]?1:0;
                }
            }
        }
        return 0;
    }
    inline std::string GetValue(std::string name, int i, int k)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                if(m_IsAdditive[j])
                {
                    if(k < m_ValuesAdditive[j].size())
                    {
                        return m_ValuesAdditive[j][k][i];
                    }else
                    {
                        return "";
                    }
                }else
                {
                    return GetValue(name,i);
                }
            }
        }
        return "";
    }
    inline double GetValueD(std::string name, int i,int k)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                if(m_IsAdditive[j])
                {
                    if(k < m_ValuesAdditive[j].size())
                    {
                        try
                        {
                            return std::stod(m_ValuesAdditive[j][k][i]);
                        }catch(...)
                        {
                            return 0.0;
                        }
                    }else
                    {
                        return 0.0;
                    }
                }else
                {
                    return GetValueD(name,i);
                }
            }
        }
        return 0;
    }

    //regular functions
    inline std::string GetValue(std::string name, int i)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                return m_Values[j][i];
            }
        }
        return "";
    }

    inline double GetValueD(std::string name, int i)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                try
                {
                    return std::stod(m_Values[j][i]);
                }catch(...)
                {
                    return 0.0;
                }

            }
        }
        return 0;
    }

    inline int GetValueI(std::string name, int i)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                return std::stod(m_Values[j][i]);
            }
        }
        return 0;
    }

    inline int GetValueB(std::string name, int i)
    {
        for(int j = 0; j < m_Names.size(); j++)
        {
            if(m_Names[j] == name)
            {
                if(is_number(m_Values[j][i]))
                {
                    return std::stod(m_Values[j][i]) > 0.5;
                }
                return (m_Values[j][i] == "True" || m_Values[j][i] == "true");
            }
        }
        return false;
    }

    inline bool ValidateS(std::string name, int i, const std::vector<std::string> &options)
    {
        std::string value = GetValue(name,i);
        for(int i = 0; i < options.size(); i++)
        {
            if(options[i] == value)
            {
                return true;
            }
        }
        return false;
    }

    inline bool ValidateIR(std::string name, int i, int value_range_min, int value_range_max)
    {
        int value = GetValueI(name,i);
        if(value >= value_range_min && value <= value_range_max)
        {
            return true;
        }
        return false;
    }

};

#endif
