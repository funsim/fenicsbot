import twitter
from time import sleep, time
from secrets import secret_dict
from parser import parser, excise

api = twitter.Api(consumer_key=secret_dict["consumer_key"],
                  consumer_secret=secret_dict["consumer_secret"],
                  access_token_key=secret_dict["access_token_key"],
                  access_token_secret=secret_dict["access_token_secret"]
)

SLEEP_TIME = 10
actually_tweet_back = True
reply_to_super_old_tweets = False

if reply_to_super_old_tweets:
    program_start_time = -1
else:
    program_start_time = time()


def tweet_image(img_fn, tweet):
    """
    Tweets an image as a reply tweet.

    :param img_fn: Filename of image to tweet.
    :param tweet: Twitter status object of tweet to reply to.
    """

    if actually_tweet_back:
        print "-----------Tweeting a solution! ------------"
        tweet_text = "@{}: I solved your problem ({})!".format(tweet.user.screen_name, excise(tweet.text).strip())[:140]
        api.PostMedia(solution_tweet_text, img_fn, 
                      in_reply_to_status_id=str(tweet.id))

def tweet_error(tweet):
    """
    Tweets an error message in reply to a tweet which did not parse.
    
    :param tweet: Tweet to reply to.
    """
    
    if actually_tweet_back:
        print ".......FEniCSbot could not come to the rescue..........."    
        error_tweet = "@{}: I failed to solve your problem...".format(tweet.user.screen_name)[:140]
        api.PostUpdate(error_tweet, in_reply_to_status_id=str(tweet.id))


last_check_id = 1
while True:
    print "--------Scanning for new tweets.-------------"
    new_mentions = api.GetSearch(term="@fenicsbot", 
                                 since_id=last_check_id)

    if len(new_mentions) == 0:
        print "--------Scan for new tweets complete.--------"
        sleep(SLEEP_TIME)
        continue

    if last_check_id == 1:
        # in first loop iteration, don't grab tweets from the beginning of time
        new_mentions = filter(lambda t: (t.created_at_in_seconds > 
                                         program_start_time), new_mentions)
    last_check_id = new_mentions[0].id


    for tweet in new_mentions:
        try:
            img_fn = parser(tweet.text)
            tweet_image(img_fn, tweet)
            
        except Exception as e :
            tweet_error(tweet)
            print "Couldn't parse: {}\n (error tweeted)".format(tweet.text)

            
    print "--------Scan for new tweets complete.--------"
    sleep(SLEEP_TIME)

