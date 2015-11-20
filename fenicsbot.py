import twitter
from time import sleep, time
from secrets import secret_dict
from parser import parser, excise

api = twitter.Api(consumer_key=secret_dict["consumer_key"],
                  consumer_secret=secret_dict["consumer_secret"],
                  access_token_key=secret_dict["access_token_key"],
                  access_token_secret=secret_dict["access_token_secret"]
)

actually_tweet_back = True
reply_to_super_old_tweets = False

if reply_to_super_old_tweets:
    program_start_time = -1
else:
    program_start_time = time()

def tweet_image(img_fn, tweet):
    """
    Tweets an image as a reply tweet.
    """
   
    # import IPython
    # IPython.embed()
    
    if not actually_tweet_back:
        print "----I want to tweet an image! Can I, please?----"    
        # import IPython
        # IPython.embed(
        
    else:
        print "-----------FEniCSbot to the rescue! ------------"    
        api.PostMedia("@{}: I solved your problem ({})!".format(tweet.user.screen_name, excise(tweet.text).strip())[:140], 
                      img_fn, in_reply_to_status_id=str(tweet.id))

def tweet_error(tweet_text, exception_text):
    """
    Tweets an error message as a friendly tweet.
    """
   
    # import IPython
    # IPython.embed()
    
    if not actually_tweet_back:
        print "----I want to tweet an error! Can I, please?----"    
        # import IPython
        # IPython.embed(
        
    else:
        print ".......FEniCSbot could not come to the rescue..........."    
        error_tweet = "@{}: I failed to solve your problem...".format(tweet.user.screen_name)
        if len(error_tweet) > 140:
            error_tweet = error_tweet[:132] + "(snip)"
        api.PostUpdate(error_tweet, in_reply_to_status_id=str(tweet.id))



last_check_id = 1
SLEEP_TIME = 10
while True:
    print "--------Scanning for new tweets.-------------"
    # print_rate_limit()
    new_mentions = api.GetSearch(term="@fenicsbot", 
                                 since_id=last_check_id)

    if len(new_mentions) == 0:
        print "--------Scan for new tweets complete.--------"
        sleep(SLEEP_TIME)
        continue

    last_check_id = new_mentions[0].id
    new_mentions = filter(lambda t: t.created_at_in_seconds > program_start_time, new_mentions)

    # print "IDs of new tweets:", map(lambda t: t.id, new_mentions)
    # print last_check_idex
    for tweet in new_mentions:
        try:
            # print tweet.text
            img_fn = parser(tweet.text)
            tweet_image(img_fn, tweet)
            
        except Exception as e :
            tweet_error(tweet.text, 
                        "{}: {}".format(e.__class__.__name__, e))
            print "Couldn't parse tweet: {}\n (error tweeted)".format(tweet.text)

            
    print "--------Scan for new tweets complete.--------"
    sleep(SLEEP_TIME)




def print_rate_limit():
    pass
    # rate_limit_status = api.GetRateLimitStatus()["resources"]["statuses"]["/statuses/mentions_timeline"]
    # print rate_limit_status
    # print "Time until rate limit reset: {}".format(rate_limit_status["reset_time"])
    # print "Number of hits remaining before reset: {}".format(rate_limit_status["remaining_hits"])
    # print "Number of hits allowed per hour: {}".format(rate_limit_status["hourly_limit"])
    # print "----------------------------------------"

    # new_mentions = api.GetMentions(since_id=last_check_id)
