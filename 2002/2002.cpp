#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstring>

#include <map>
#include <string>

using namespace std;

struct TUser
{
    string _password;
    bool _loggedIn;

    TUser()
        : _loggedIn(false)
    {

    }
};

typedef map<string, TUser> TUsers;

int main()
{
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
#endif

    int n;
    scanf("%d", &n);
    TUsers users;
    for (int i = 0; i < n; ++i)
    {
        char sCommand[1000];
        scanf("%s", sCommand);
        if (!strcmp(sCommand, "register"))
        {
            char username[1000];
            char password[1000];
            scanf("%s%s", username, password);
            string sUsername(username);
            if (users.find(sUsername) != users.end())
            {
                printf("fail: user already exists\n");
            }
            else
            {
                TUser user;
                user._password = password;
                users[sUsername] = user;
                printf("success: new user added\n");
            }
        }
        else if (!strcmp(sCommand, "login"))
        {
            char username[1000];
            char password[1000];
            scanf("%s%s", username, password);
            string sUsername(username);
            if (users.find(sUsername) != users.end())
            {
                TUser& user = users[sUsername];
                if (user._password == password)
                {
                    if (user._loggedIn)
                    {
                        printf("fail: already logged in\n");
                    }
                    else
                    {
                        user._loggedIn = true;
                        printf("success: user logged in\n");
                    }
                }
                else
                {
                    printf("fail: incorrect password\n");
                }
            }
            else
            {
                printf("fail: no such user\n");
            }
        }
        else if (!strcmp(sCommand, "logout"))
        {
            char username[1000];
            scanf("%s", username);
            string sUsername(username);
            if (users.find(sUsername) != users.end())
            {
                TUser& user = users[sUsername];
                if (!user._loggedIn)
                {
                    printf("fail: already logged out\n");
                }
                else
                {
                    user._loggedIn = false;
                    printf("success: user logged out\n");
                }
            }
            else
            {
                printf("fail: no such user\n");
            }
        }
    }
    return 0;
}