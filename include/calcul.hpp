//////////////////////////////////////////////////////////////
//author lovejava (alias lovecpp, real name : Yann Hamdaoui)//
//////////////////////////////////////////////////////////////
//This source code is under the GNU/GPL license.            //
//////////////////////////////////////////////////////////////

#ifndef _CALCUL_HPP_
#define _CALCUL_HPP_

#include <iostream>
#include <map>
#include <boost/algorithm/string.hpp>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

/**********************************************
 *d�finit une fonctions (et ses sous fonctions)
 *capable d'�valuer une expression math�matique sous forme de cha�ne
 *et de renvoyer le r�sultat
 */

namespace calcul
{
    /**
     *exception signalant que la syntaxe de l'expression pass�e en param�tre
     *n'est pas valide
     */
    class syntax_exception : public exception
    {
        protected:
            string message;
            int position; //position du caract�re ayant provoqu� l'erreur
            //expression math�matique compl�te
            string chaine;

        public:
            syntax_exception() throw() { };
            syntax_exception(const string& message, const string& chaine = "", int position = string::npos) throw();
            syntax_exception(const syntax_exception& source) throw();

            const char* what() const throw();

            syntax_exception& operator=(const syntax_exception& source) throw();

            virtual ~syntax_exception() throw() { };
    };

    /**
     *exception lev�e lors d'une division par z�ro
     */
    class division_by_zero_exception : public exception
    {
        public:
            division_by_zero_exception() throw() { };
            division_by_zero_exception(const division_by_zero_exception& source) throw() { };

            const char* what() const throw();

            division_by_zero_exception& operator=(const division_by_zero_exception& source) throw() { };

            virtual ~division_by_zero_exception() throw() { };
    };

    /**
     *exception lev�e lorsqu'une variable ou une fonction utilis�e dans un calcul n'est pas d�finie
     */
    class undefined_symbol_exception : public exception
    {
        protected:
            string symbole; //le nom de la variable/fonction ind�finie

        public:
            undefined_symbol_exception(const string& symbole = "?") throw();
            undefined_symbol_exception(const undefined_symbol_exception& source) throw();

            const char* what() const throw();

            undefined_symbol_exception& operator=(const undefined_symbol_exception& source) throw();

            virtual ~undefined_symbol_exception() throw() { };
    };

    /**
     *Classe qui permet d'envelopper les exceptions qui ont sont lev�es lorsque l'�valuateur
     *parse une expression de fonction. Par exemple, si l'utilisateur d�finit la fonction f : (x^2 + 2x + 5)/(x-1)
     *et qu'il tape 2*_f(1)/2, il peut ne pas comprendre d'o� vient l'erreur si le seul message qu'il obtient est "division by zero".
     *C'est encore plus un probl�me lorsqu'une erreur de syntaxe est pr�sente dans la d�claration de la fonction.
     *Cette exception permet donc de palier � ce probl�mes en encapsulant l'erreur de base et en rajoutant dans le message d'erreur final
     *la provenance de l'erreur : par exemple [in function f]division by zero
     */
    class fonction_error_wrapper_exception : public exception
    {
        protected:
            string nom; //le nom de la fonction d'o� provient l'erreur
            string message_origine; //le message d'erreur d'origine

        public:
            fonction_error_wrapper_exception(const exception& erreur, const string& nom = "?") throw();
            fonction_error_wrapper_exception(const fonction_error_wrapper_exception& source) throw();

            const char* what() const throw();

            fonction_error_wrapper_exception& operator=(const fonction_error_wrapper_exception& source) throw();

            virtual ~fonction_error_wrapper_exception() throw() { };
    };

    #define ALGORITHMIC_FUNCTION number_type (*)(number_type)

    /**
     *Classe qui �value une expression pour la parser sous forme de nombre.
     *le param�tre template number_type permet de l'utiliser aussi bien avec
     *des int, des doubles, des float, etc.
     *le constructeur seul ne fait qu'initialiser l'objet, le parsage de l'expression
     *se fait lors de la conversion en objet de type number_type (voir operator number_type())
     *ou de l'appel explicite � eval()
     */
    template<typename number_type, typename algorithmic_function = number_type (*)(number_type)>
    class expression_evaluator
    {
        /*public:
            typedef (number_type (*)(number_type)) algorithmic_function;*/

        protected:
            string expression; //l'expression � parser
            unsigned int index; //index utilis� lors du parsage
            map<char, number_type> variables; //variables pouvant �tre utilis�es pour des expressions
            map<string, expression_evaluator<number_type> > fonctions; //fonctions d�finies par l'utilisateur
            map<string, algorithmic_function> fonctions_algo; //fonctions algorithmiques et pr�d�finies

            //m�thodes d'analyse et de parcours de la cha�ne (voir le corps des m�thodes
            //pour des commentaires plus d�taill�s)
            number_type somme(); //+ et -
            number_type prior(); //*, / et %
            number_type puissance(); // ^
            number_type terme(); //analyse les nombres � proprement dit

            bool more(); //si l'on est arriv� � la fin de la cha�ne
            /*Fonction qui renvoie true s'il y a un signe multipli� implicite � la position courante,
             *comme par exemple dans 2x + 3ab(2 + 3z)d, c'est � dire lorsque des lettres, des nombres et/ou des parenth�ses
             *se suivent sans op�rateur*/
            bool implicit_mult();

        public:
            expression_evaluator() { };
            expression_evaluator(const string& expression);
            expression_evaluator(const expression_evaluator& source);

            string get_expression() const;
            /*cette fonction renvoie la valeur d'une variable, ou 0 si elle n'existe pas
             *pour v�rifier qu'une variable existe, utiliser is_var_defined() au pr�alable, car une variable
             *existante peut tr�s bien avoir pour valeur 0
             */
            number_type get_var(char variable) const;
            /*
             *renvoie l'expression d'une fonction d�finie par l'utilisateur
             */
            string get_function(const string& nom) const;
            /*
             *renvoie un pointeur de fonction vers une fonction "algorithmique" d�finie par l'utilisateur
             *ou pr�d�finie (sinus, cosinus, etc...). Par algorithmique, j'entend une fonction qui ne peut �tre
             *d�finie seulement par une expression math�matique comme '2x+3', par exemple la fonction cosinus, partie enti�re
             *ou encore racine carr�e.
             */
            algorithmic_function get_algorithmic_function(const string& nom) const;

            void set_expression(const string& expression);
            //permet d'assigner une valeur � une variable
            void set_var(char variable, number_type valeur);
            //permet d'assigner � une fonction son expression
            void set_function(const string& nom, const string& expression);
            //permet d'assigner � une fonction son algorithme
            void set_algorithmic_function(const string& nom, algorithmic_function fonction);

            bool is_var_defined(char variable) const;
            bool is_function_defined(const string& nom) const;

            //�value l'expression et renvoie un nombre
            number_type eval();
            //v�rifie la syntaxe. L�ve une syntax_exception si celle-ci n'est pas valide
            void check_syntax() const throw(syntax_exception);

            //m�thode de classification des caract�res
            static bool is_digit(char caractere);
            static bool is_alpha(char caractere);
            static bool is_alphanum(char caractere);
            static bool is_operator(char caractere);
            static bool is_forbidden(char caractere);

            expression_evaluator& operator=(const expression_evaluator& source);
            operator number_type();

            ~expression_evaluator();
    };

    //fonctions d'adaptation permettant d'utiliser les fonctions usuelles comme cosinus, log, etc...
    //avec n'importe quel type g�n�rique number_type en effectuant les conversions n�cessaires :
    template<typename number_type>
    number_type _abs_adapter(number_type parametre)
    {
        return (number_type)fabs((double)parametre);
    }


    template<typename number_type>
    number_type _acos_adapter(number_type parametre)
    {
        return (number_type)acos((double)parametre);
    }

    template<typename number_type>
    number_type _asin_adapter(number_type parametre)
    {
        return (number_type)asin((double)parametre);
    }

    template<typename number_type>
    number_type _atan_adapter(number_type parametre)
    {
        return (number_type)atan((double)parametre);
    }


    template<typename number_type>
    number_type _ceil_adapter(number_type parametre)
    {
        return (number_type)ceil((double)parametre);
    }


    template<typename number_type>
    number_type _exp_adapter(number_type parametre)
    {
        return (number_type)exp((double)parametre);
    }


    template<typename number_type>
    number_type _floor_adapter(number_type parametre)
    {
        return (number_type)floor((double)parametre);
    }


    template<typename number_type>
    number_type _log_adapter(number_type parametre)
    {
        return (number_type)log((double)parametre);
    }

    template<typename number_type>
    number_type _log10_adapter(number_type parametre)
    {
        return (number_type)log10((double)parametre);
    }


    template<typename number_type>
    number_type _cos_adapter(number_type parametre)
    {
        return (number_type)cos((double)parametre);
    }

    template<typename number_type>
    number_type _cosh_adapter(number_type parametre)
    {
        return (number_type)cosh((double)parametre);
    }

    template<typename number_type>
    number_type _sin_adapter(number_type parametre)
    {
        return (number_type)sin((double)parametre);
    }

    template<typename number_type>
    number_type _sinh_adapter(number_type parametre)
    {
        return (number_type)sinh((double)parametre);
    }

    template<typename number_type>
    number_type _tan_adapter(number_type parametre)
    {
        return (number_type)tan((double)parametre);
    }

    template<typename number_type>
    number_type _tanh_adapter(number_type parametre)
    {
        return (number_type)tanh((double)parametre);
    }

    ///////////////////////////////////////////////////////

    syntax_exception::syntax_exception(const string& message, const string& chaine, int position) throw() : exception()
    {
        this->message = message;
        this->chaine = chaine;
        this->position = position;
    }

    syntax_exception::syntax_exception(const syntax_exception& source) throw() : exception()
    {
        *this = source;
    }

    const char* syntax_exception::what() const throw()
    {
        string temp(message);
        //si un num�ro de caract�re est sp�cifi� pour l'erreur syntaxique, on affiche le caract�re en question
        //et sa position dans l'expression apr�s le message d'erreur
        if(position != string::npos)
        {
            ostringstream sortie;
            sortie << position+1;
            temp.append("(caract�re \'").append(chaine, position, 1).append("\', en position ").append(sortie.str()).append(")");
        }

        return temp.c_str();
    }

    syntax_exception& syntax_exception::operator=(const syntax_exception& source) throw()
    {
        message = source.message;
        chaine = source.chaine;
        position = source.position;
    }

    ///////////////////////////////////////////////////////////////////////////

    const char* division_by_zero_exception::what() const throw()
    {
        return "division par z�ro";
    }

    ///////////////////////////////////////////////////////////////////////////

    undefined_symbol_exception::undefined_symbol_exception(const string& symbole) throw()
    {
        this->symbole = symbole;
    }

    undefined_symbol_exception::undefined_symbol_exception(const undefined_symbol_exception& source) throw()
    {
        *this = source;
    }

    const char* undefined_symbol_exception::what() const throw()
    {
        return string("\'").append(symbole).append("\' n'est pas d�fini.").c_str();
    }

    undefined_symbol_exception& undefined_symbol_exception::operator=(const undefined_symbol_exception& source) throw()
    {
        symbole = source.symbole;
    }

    ///////////////////////////////////////////////////////////////////////////

    fonction_error_wrapper_exception::fonction_error_wrapper_exception(const exception& erreur, const string& nom) throw()
    {
        this->message_origine = erreur.what();
        this->nom = nom;
    }

    fonction_error_wrapper_exception::fonction_error_wrapper_exception(const fonction_error_wrapper_exception& source) throw()
    {
        *this = source;
    }

    const char* fonction_error_wrapper_exception::what() const throw()
    {
        return string("[in function ").append(nom).append("]").append(message_origine).c_str();
    }

    fonction_error_wrapper_exception& fonction_error_wrapper_exception::operator=(const fonction_error_wrapper_exception& source) throw()
    {
        nom = source.nom;
        message_origine = source.message_origine;
    }

    ///////////////////////////////////////////////////////////////////////////

    template<typename number_type, typename algorithmic_function>
    expression_evaluator<number_type, algorithmic_function>::expression_evaluator(const string& expression)
    {
        this->expression = expression;

        //ajout des fonctions math�matiques de base
        fonctions_algo["abs"] = &_abs_adapter<number_type>;

        fonctions_algo["acos"] = &_acos_adapter<number_type>;
        fonctions_algo["asin"] = &_asin_adapter<number_type>;
        fonctions_algo["atan"] = &_atan_adapter<number_type>;

        fonctions_algo["ceil"] = &_ceil_adapter<number_type>;

        fonctions_algo["exp"] = &_exp_adapter<number_type>;

        fonctions_algo["floor"] = &_floor_adapter<number_type>;

        fonctions_algo["log"] = &_log_adapter<number_type>;
        fonctions_algo["log10"] = &_log10_adapter<number_type>;

        fonctions_algo["cos"] = &_cos_adapter<number_type>;
        fonctions_algo["cosh"] = &_cosh_adapter<number_type>;
        fonctions_algo["sin"] = &_sin_adapter<number_type>;
        fonctions_algo["sinh"] = &_sinh_adapter<number_type>;
        fonctions_algo["tan"] = &_tan_adapter<number_type>;
        fonctions_algo["tanh"] = &_tanh_adapter<number_type>;
    }

    template<typename number_type, typename algorithmic_function>
    expression_evaluator<number_type, algorithmic_function>::expression_evaluator(const expression_evaluator<number_type, algorithmic_function>& source)
    {
        *this = source;
    }

    template<typename number_type, typename algorithmic_function>
    number_type expression_evaluator<number_type, algorithmic_function>::somme()
    {
        /*on commence par s�parer les termes via les signes + et -,
         *ceux qui ont la plus faible priorit�. On calculera d'abord
         *ces termes ind�pendamment pour ensuite leur appliquer, en dernier lieu,
         *l'op�rateur + ou -.
         */


        //� l'appel de la fonction, on consid�re que le premier �l�ment est forc�ment
        //soit un nombre, soit une parenth�se : on appelle donc la m�thode terme()
        //via la m�thode prior() pour r�cup�rer ce premier nombre.
        number_type retour = prior();

        //une fois obtenu un nombre, on cherche un signe + ou -
        //si on trouve, on appel les m�thodes sous-jacente pour calculer
        //l'ensemble de la valeur qui sera ajout�e ou retranch�e au final
        while(more() && (expression[index] == '+' || expression[index] == '-'))
        {
            ++index; //on avance d'un caract�re
            if(expression[index-1] == '+')
                retour += prior();
            else
                retour -= prior();
        }

        /*si on arrive � la fin d'une parenth�se, on sort d'une sous-expression.
         *On avance donc la position courante pour que, lorsque l'on rend la main �
         *la fonction terme(), l'analyse de l'expression ne s'arr�te pas pr�matur�ment
         *� cause du fait que comme la parenth�se n'est ni un op�rateur, ni un nombre,
         *les fonctions rendraient la main une � une sans effectuer aucun traitement suppl�mentaire
         *Lorsqu'on enl�ve cette ligne par exemple, (5+(5+5)*2)*2 font 15 (en fait, seulement les
         *trois premiers 5 on �t� analys�s et ajout�s) et non pas 50.
         */
        if(more() && expression[index] == ')')
            ++index;

        return retour;
    }

    template<typename number_type, typename algorithmic_function>
    number_type expression_evaluator<number_type, algorithmic_function>::prior()
    {
        number_type retour = puissance(); //on obtient le premier nombre de l'expression via puissance()

        /*une fois obtenu un nombre, on cherche un signe *, / ou %
         *on continue jusqu'� que le signe suivant ne soit ni un nombre ni un * ni un / ni un %
         *pour �tre s�r d'avoir effectu� toutes les divisions/multiplications de l'expression
         *ou sous-expression. Ainsi la priorit� est respect�e et on rend la main � la fonction somme
         */
        bool implicite = false; //on doit savoir si la multiplication �tait implicite ou non pour avancer ou pas index
        while( more() && (expression[index] == '*' || expression[index] == '/' || expression[index] == '%' || (implicite = implicit_mult())) )
        {
            if(!implicite)
                ++index;
            if(expression[index-1] == '/')
            {
                number_type temp = puissance();
                if(temp == 0)
                    throw division_by_zero_exception();

                retour /= temp;
            }
            else if(expression[index-1] == '%')
                retour = (number_type) ((long long int)retour % (long long int)puissance());
            else
                retour *= puissance();

            /*on r�initialise implicite pour les prochains passages, car si expression[index] vaut '*' par exemple, l'instruction
             *implicite = implicit_mult() ne sera pas ex�cut�e ; donc si au pr�c�dent passage implicite valait true, il garderait cette valeur
             *si on ne le remettait pas � false
             */
            implicite = false;
        }

        return retour;
    }

    template<typename number_type, typename algorithmic_function>
    number_type expression_evaluator<number_type, algorithmic_function>::puissance()
    {
        /*fonction qui regarde si dans l'expression ou sous-expression en cours
         *il n'y a pas de signe '^'. Elle est apr�s la fonction prior() dans le cha�nage
         *des appels car l'op�rateur d'�l�vation � la puissance est le plus prioritaire
         */

        number_type retour = terme();

        while(more() && (expression[index] == '^'))
        {
            ++index;
            retour = (number_type)pow((double)retour, (double)terme());
        }

        return retour;
    }

    template<typename number_type, typename algorithmic_function>
    number_type expression_evaluator<number_type, algorithmic_function>::terme()
    {
        //c'est la fonction qui lit r�ellement les nombres
        //si on tombe sur une parenth�se ouvrante, on r�-analyse
        //cette sous-expression comme une nouvelle expression apparenti�re
        if(expression[index] == '(')
        {
            ++index;
            return somme();
        }
        else
        {
            number_type nombre;
            unsigned int index_depart = index;
            bool moins_unaire = false; //s'il y avait un moins unaire au d�but de l'expression

            //gestion du moins unaire : si l'on rencontre un moins lorsque les fonctions
            //appellantes veulent un nombre, on consid�re que c'est un moins unaire
            //on met un bool�en � true puis on avance l'index courant
            if(expression[index] == '-')
            {
                moins_unaire = true;
                ++index;
                ++index_depart;
            }

            if(expression_evaluator<number_type, algorithmic_function>::is_operator(expression[index]))
                throw syntax_exception("plusieurs op�rateurs se suivent", expression, index);

            if(expression[index] == '_')
            {
                //on avance la position courante jusqu'� la parenth�se ouvrante
                ++index;
                ++index_depart;

                for(; more() && expression[index] != '('; ++index) ;
                //on doit r�cup�rer la valeur d'une fonction
                if(!more())
                    throw syntax_exception("Un nom de fonction doit �tre suivi d'une parenth�se ouvrante", expression, index_depart);

                string nom_fonction = expression.substr(index_depart, index-index_depart);
                typename map<string, /*expression_evaluator<number_type, algorithmic_function>::algorithmic_function*/ALGORITHMIC_FUNCTION>::iterator it1 = fonctions_algo.find(nom_fonction);
                typename map<string, expression_evaluator<number_type, algorithmic_function> >::iterator it2 = fonctions.find(nom_fonction);

                ++index; //on passe apr�s la parenth�se
                number_type x_value = somme(); //on �value ind�pendamment l'expression qui suit la fonction (c'est � dire la valeur de x)

                if(it1 != fonctions_algo.end())
                {
                    try {
                        nombre = (*(*it1).second)(x_value);
                    }
                    //il faut convertir les messages d'erreurs pour bien sp�cifier que l'erreur se trouve dans la d�finition de la fonction,
                    //et non pas dans la cha�ne analys�e en ce moment
                    catch(const exception& e)
                    {
                        throw fonction_error_wrapper_exception(e, nom_fonction);
                    }
                }
                else if(it2 != fonctions.end())
                {

                    //on met � jour les d�finitions de variables et de fonctions pour que la suite tout �a etc.
                    fonctions[nom_fonction].variables = variables;
                    fonctions[nom_fonction].fonctions = fonctions;
                    fonctions[nom_fonction].set_var('x', x_value);

                    try {
                        nombre = fonctions[nom_fonction].eval();
                    }
                    //il faut convertir les messages d'erreurs pour bien sp�cifier que l'erreur se trouve dans la d�finition de la fonction,
                    //et non pas dans la cha�ne analys�e en ce moment
                    catch(const exception& e)
                    {
                        throw fonction_error_wrapper_exception(e, nom_fonction);
                    }
                }
                else
                    throw undefined_symbol_exception(nom_fonction);
            }
            else if(expression_evaluator<number_type, algorithmic_function>::is_alpha(expression[index]))
            {
                char variable = expression[index];
                if(!is_var_defined(variable))
                    throw undefined_symbol_exception(string(&variable, 1));

                number_type valeur = get_var(variable);
                nombre = valeur;
                ++index;
            }
            else
            {
                //sinon, on lit le nombre � la position courante
                for(; more() && is_digit(expression[index]); ++index) ;

                string nombre_str(expression, index_depart, index-index_depart);
                istringstream entree(nombre_str);
                entree >> nombre;
            }

            if(moins_unaire)
                return -nombre;
            else
                return nombre;
        }
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::more()
    {
        return index < expression.length();
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::implicit_mult()
    {
        /*V�rifie s'il n'y a pas de signe multiplier implicite
         *avec des variables (comme par exemple 2a(3+2) ou encore (x+2)(x-3)bc
         *v�rifications : si un chiffre et une lettre se suivent (ou l'inverse)
         *                si deux lettres se suivent
         *                si un chiffre/une lettre suivent une parenth�se fermante (ou pr�c�dent une parenth�se ouvrante)
         *                si une parenth�se fermante et une parenth�se ouvrante se suivent
         */
        return ( (expression_evaluator<number_type, algorithmic_function>::is_digit(expression[index-1]) && expression_evaluator<number_type, algorithmic_function>::is_alpha(expression[index])) ||
                 (expression_evaluator<number_type, algorithmic_function>::is_alpha(expression[index-1]) && expression_evaluator<number_type, algorithmic_function>::is_digit(expression[index])) ||
                 (expression_evaluator<number_type, algorithmic_function>::is_alpha(expression[index-1]) && expression_evaluator<number_type, algorithmic_function>::is_alpha(expression[index])) ||
                 (expression_evaluator<number_type, algorithmic_function>::is_alphanum(expression[index-1]) && expression[index] == '(') ||
                 (expression[index-1] == ')' && expression_evaluator<number_type, algorithmic_function>::is_alphanum(expression[index])) ||
                 (expression[index-1] == ')' && expression[index] == '(') );
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::is_digit(char caractere)
    {
        return (caractere == '.' ||
                (caractere >= '0' && caractere <= '9'));
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::is_alpha(char caractere)
    {
        return (caractere >= 'a' && caractere <= 'z') || (caractere >= 'A' && caractere <= 'Z') || caractere == '_';
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::is_alphanum(char caractere)
    {
        return is_digit(caractere) || is_alpha(caractere);
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::is_operator(char caractere)
    {
        return (caractere == '+' || caractere == '-' || caractere == '%'
                || caractere == '*' || caractere == '/' || caractere == '^');
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::is_forbidden(char caractere)
    {
        return (!is_digit(caractere) && !is_operator(caractere) && !is_alpha(caractere)
                && caractere != '(' && caractere != ')' && caractere != ' ');
    }

    template<typename number_type, typename algorithmic_function>
    string expression_evaluator<number_type, algorithmic_function>::get_expression() const
    {
        return expression;
    }

    template<typename number_type, typename algorithmic_function>
    number_type expression_evaluator<number_type, algorithmic_function>::get_var(char variable) const
    {
        typename map<char, number_type>::const_iterator it = variables.find(variable);

        if(it == variables.end())
            return 0;
        else
            return (*it).second;
    }

    template<typename number_type, typename algorithmic_function>
    string expression_evaluator<number_type, algorithmic_function>::get_function(const string& nom) const
    {
        typename map<string, expression_evaluator<number_type, algorithmic_function> >::const_iterator it = fonctions.find(nom);

        if(it == fonctions.end())
            return "";
        else
            return (*it).second.get_expression();
    }

    template<typename number_type, typename algorithmic_function>
    algorithmic_function expression_evaluator<number_type, algorithmic_function>::get_algorithmic_function(const string& nom) const
    {
        typename map<string, /*expression_evaluator<number_type, algorithmic_function>::algorithmic_function*/ALGORITHMIC_FUNCTION>::const_iterator it = fonctions.find(nom);

        if(it == fonctions.end())
            return NULL;
        else
            return (*it).second;
    }

    template<typename number_type, typename algorithmic_function>
    void expression_evaluator<number_type, algorithmic_function>::set_expression(const string& expression)
    {
        this->expression = expression;
    }

    template<typename number_type, typename algorithmic_function>
    void expression_evaluator<number_type, algorithmic_function>::set_var(char variable, number_type valeur)
    {
        variables[variable] = valeur;
    }

    template<typename number_type, typename algorithmic_function>
    void expression_evaluator<number_type, algorithmic_function>::set_function(const string& nom, const string& expression)
    {
        //on cr�er une nouvelle fonction seulement si elle n'est pas d�j� d�finie parmis les fonctions algorithmiques ou pr�d�finies
        //pour �viter les conflits
        if(fonctions_algo.find(nom) == fonctions_algo.end())
            fonctions[nom] = expression;
    }

    template<typename number_type, typename algorithmic_function>
    void expression_evaluator<number_type, algorithmic_function>::set_algorithmic_function(const string& nom, algorithmic_function fonction)
    {
        //si la fonction existe d�j� parmis les fonctions d�finies par l'utilisateur, on la supprime de la liste avant de l'ins�rer
        //parmis les fonctions algorithmiques
        typename map<string, expression_evaluator<number_type, algorithmic_function> >::iterator it = fonctions.find(nom);
        if(it != fonctions.end())
            fonctions.erase(it);
        fonctions_algo[nom] = fonction;
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::is_var_defined(char variable) const
    {
        return variables.find(variable) != variables.end();
    }

    template<typename number_type, typename algorithmic_function>
    bool expression_evaluator<number_type, algorithmic_function>::is_function_defined(const string& nom) const
    {
        return fonctions.find(nom) != fonctions.end() || fonctions_algo.find(nom) != fonctions_algo.end();
    }

    template<typename number_type, typename algorithmic_function>
    number_type expression_evaluator<number_type, algorithmic_function>::eval()
    {
        if(expression.find(' ') != string::npos)
            boost::algorithm::erase_all(expression, " "); //suppression des espaces

        check_syntax();
        index = 0;
        return somme();
    }

    template<typename number_type, typename algorithmic_function>
    void expression_evaluator<number_type, algorithmic_function>::check_syntax() const throw(syntax_exception)
    {
        if((expression_evaluator<number_type, algorithmic_function>::is_operator(expression[0]) && expression[0] != '-')
           || expression_evaluator<number_type, algorithmic_function>::is_operator(expression[expression.length()-1]))
            throw syntax_exception("l'expression commence ou se termine par un op�rateur");

        register int i;
        register int compteur = 0;

        for(i = 0; i < expression.length(); ++i)
        {
            if(expression_evaluator<number_type, algorithmic_function>::is_forbidden(expression[i]))
                throw syntax_exception("caract�re interdit", expression, i);
            if(expression[i] == '(')
                ++compteur;
            if(expression[i] == ')')
                --compteur;
        }

        if(compteur != 0)
            throw syntax_exception("le nombre de parenth�ses ouvrantes ne correspond pas au nombres de parenth�ses fermantes.");
    }

    template<typename number_type, typename algorithmic_function>
    expression_evaluator<number_type, algorithmic_function>& expression_evaluator<number_type, algorithmic_function>::operator=(const expression_evaluator<number_type, algorithmic_function>& source)
    {
        expression = source.expression;
        variables = source.variables;
        fonctions = source.fonctions;
        fonctions_algo = source.fonctions_algo;
    }

    template<typename number_type, typename algorithmic_function>
    expression_evaluator<number_type, algorithmic_function>::operator number_type()
    {
        return eval();
    }

    template<typename number_type, typename algorithmic_function>
    expression_evaluator<number_type, algorithmic_function>::~expression_evaluator()
    {
    }
}

#endif /* ndef _CALCUL_HPP_ */
