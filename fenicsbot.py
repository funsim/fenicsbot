import twitter
from time import sleep, time
from secrets import secret_dict
from parser import parser

api = twitter.Api(consumer_key=secret_dict["consumer_key"],
                  consumer_secret=secret_dict["consumer_secret"],
                  access_token_key=secret_dict["access_token_key"],
                  access_token_secret=secret_dict["access_token_secret"]
)

def tweet_image(img_fn, reply_to_id):
    """
    Tweets an image as a reply tweet.
    """
    print "----I want to tweet an image! Can I, please?----"    
    # api.PostMedia("I solved your problem!", img_fn, 
    #               in_reply_to_status_id=reply_to_id)


last_check_id = 1
program_start_time = -1
# program_start_time = time()
SLEEP_TIME = 10
while True:
    print "--------Scanning for new tweets.-------------"
    # print_rate_limit()
    new_mentions = api.GetSearch(term="@fenicsbot", since_id=last_check_id)

    if len(new_mentions) == 0:
        sleep(SLEEP_TIME)
        print "--------Scan for new tweets complete.--------"
        continue

    last_check_id = new_mentions[0].id
    new_mentions = filter(lambda t: t.created_at_in_seconds > program_start_time, new_mentions)

    # print "IDs of new tweets:", map(lambda t: t.id, new_mentions)
    # print last_check_idex
    for tweet in new_mentions:
        try:
            # print tweet.text
            img_fn = parser(tweet.text)
            tweet_image(img_fn, tweet.id)
            
        except Exception as e :
            print "Couldn't parse tweet: {}".format(tweet.text)
            print "{}: {}".format(e.__class__.__name__, e)
            
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
